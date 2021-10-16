import os
import sys
import glob
import numpy as np
from jet_image import JetImage, TwinJetImage
from vlbi_utils import get_transverse_profile
sys.path.insert(0, '../ve/vlbi_errors')
from uv_data import UVData, downscale_uvdata_by_freq
from spydiff import clean_difmap
from from_fits import create_clean_image_from_fits_file, create_image_from_fits_file, create_model_from_fits_file
from image_ops import spix_map
from image import plot as iplot

n_rep = 12
mhd_image = "/home/ilya/github/mhd_transfer/Release/m1s10g2b123.971372r0.000369_psi_1.000000_dpsi_0.015000_jet_image_i_15.4.txt"
save_dir = "/home/ilya/data/M87Lesha/artMOJAVE/stack/newmod"
data_dir = "/home/ilya/data/M87Lesha/artMOJAVE/stack"
uvfits_files = sorted(glob.glob(os.path.join(data_dir, "1228+126*.uvf")))

path_to_script = "/home/ilya/github/stackemall/external_scripts/last/final_clean_rms_last"
common_mapsize = (1024, 0.1)
# MOJAVE receipt using DEC
common_beam = (0.85, 0.85, 0)
noise_scale_factor = 1.0
z = 0.00436
n_along = 800
n_across = 150
lg_pixel_size_mas_min = np.log10(0.025)
lg_pixel_size_mas_max = np.log10(0.05)
jm = JetImage(z=z, n_along=n_along, n_across=n_across,
              lg_pixel_size_mas_min=lg_pixel_size_mas_min, lg_pixel_size_mas_max=lg_pixel_size_mas_max,
              jet_side=True, rot=np.deg2rad(-107.0))
jm.load_image_stokes("I", mhd_image, scale=0.5)

all_images = dict()
all_stacks = list()

for i_rep in range(n_rep):
    images = list()
    bad_epochs = list()
    for i, uvfits_file in enumerate(uvfits_files[:]):
        print("==============================================")
        print("Working with ", uvfits_file)
        print("==============================================")
        try:
            uvdata = UVData(uvfits_file)
        # Some fucking weird epochs
        except:
            bad_epochs.append(uvfits_file)
            continue
        noise = uvdata.noise(average_freq=False, use_V=False)
        uvdata.zero_data()
        # If one needs to decrease the noise this is the way to do it
        for baseline, baseline_noise_std in noise.items():
            noise.update({baseline: noise_scale_factor*baseline_noise_std})

        uvdata.substitute([jm])
        uvdata.noise_add(noise)
        uvdata.save("arty.uvf", rewrite=True, downscale_by_freq=downscale_uvdata_by_freq(uvdata))
        #
        clean_difmap(fname="arty.uvf",
                     outfname="arty_cc.fits", outpath=save_dir, stokes="i",
                     mapsize_clean=common_mapsize, path_to_script=path_to_script,
                     show_difmap_output=True,
                     beam_restore=common_beam,
                     dfm_model=os.path.join(save_dir, "arty.mdl"))
        # clean_difmap(fname=uvfits_file,
        #              outfname="arty_cc.fits", outpath=save_dir, stokes="i",
        #              mapsize_clean=common_mapsize, path_to_script=path_to_script,
        #              show_difmap_output=True,
        #              beam_restore=common_beam,
        #              dfm_model=os.path.join(save_dir, "arty.mdl"))

        image = create_image_from_fits_file(os.path.join(save_dir, "arty_cc.fits")).image
        images.append(image)
        os.unlink(os.path.join(save_dir, "arty_cc.fits"))

    # Keep all individual images of current stack realization
    all_images[i_rep] = images

    from astropy.stats import mad_std
    good_images = [im for im in images if mad_std(im) < 0.01]
    stack_image = np.sum(good_images, axis=0)/len(good_images)
    all_stacks.append(stack_image)
    get_transverse_profile(stack_image, PA=107, nslices=200, plot_zobs_min=1, plot_zobs_max=20, beam=common_beam,
                           pixsize_mas=0.1, treat_as_numpy_array=True, save_figs=True, save_dir=save_dir,
                           save_prefix="M87_MDH_{}".format(str(i_rep+1).zfill(2)))
    print("Bad epochs : ", bad_epochs)


# Plot several stacks
fig, fig_res = get_transverse_profile(all_stacks[0], PA=107, nslices=200, plot_zobs_min=1, plot_zobs_max=20, beam=common_beam,
                                      pixsize_mas=0.1, treat_as_numpy_array=True, save_figs=False, alpha=0.2)
for im in all_stacks[1:]:
    fig, fig_res = get_transverse_profile(im, PA=107, nslices=200, plot_zobs_min=1, plot_zobs_max=20, beam=common_beam,
                                 pixsize_mas=0.1, treat_as_numpy_array=True, fig=fig, fig_res=fig_res, save_figs=False, alpha=0.2)
fig.savefig(os.path.join(save_dir, "test.png"), bbox_inches="tight", dpi=300)
fig_res.savefig(os.path.join(save_dir, "testres.png"), bbox_inches="tight", dpi=300)