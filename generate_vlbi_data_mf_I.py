import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import sys
from jet_image import JetImage, TwinJetImage
from vlbi_utils import (find_image_std, find_bbox, pol_mask, correct_ppol_bias, convert_difmap_model_file_to_CCFITS,
                        rotate_difmap_model)
sys.path.insert(0, '../ve/vlbi_errors')
from uv_data import UVData
from spydiff import clean_difmap
from from_fits import create_clean_image_from_fits_file, create_image_from_fits_file
from image_ops import spix_map
from image import plot as iplot

only_make_pics = False

# Directory to save all
# Original directory
# save_dir = "/home/ilya/data/M87Lesha/art"
# save_dir = "/home/ilya/data/M87Lesha/artMOJAVE/central_ridge"
# Saving intermediate files
save_dir = "/home/ilya/data/mf/spix/2ridges"
# Saving final images
save_dir2 = "/home/ilya/data/mf/spix/2ridges"
if save_dir2 is None:
    save_dir2 = save_dir
# Observed frequencies of simulations
freqs_obs_ghz = [8.1, 15.4]
# freqs_obs_ghz = [15.4]

# S ~ nu^{+alpha}
alpha_true = -0.5

# Multiplicative factor for noise added to model visibilities.
noise_scale_factor = 0.5
# Common size of the map and pixel size (mas)
common_mapsize = (1024, 0.1)
# Common beam size (mas, mas, deg)
common_beam = (1.5, 1.5, 0)
template_x_ccimage = create_clean_image_from_fits_file("/home/ilya/data/M87Lesha/to_ilya/X_template_beam.fits")
common_beam = template_x_ccimage.beam

mhdrt_code = "m1s10g2b123.971372r0.000369_psi_1.000000_dpsi_0.015000"
# Original directory
# jetpol_files_directory = "/home/ilya/github/mhd_transfer/results/n5000"
# jetpol_files_directory = "/home/ilya/github/mhd_transfer/Release"
# jetpol_files_directory = "/home/ilya/data/M87Lesha/mhd_images_bkp/3ridges"
jetpol_files_directory = "/home/ilya/data/M87Lesha/mhd_images_bkp/2ridges"

# C++ code run parameters
z = 0.00436
n_along = 800
n_across = 150
lg_pixel_size_mas_min = np.log10(0.025)
lg_pixel_size_mas_max = np.log10(0.05)

path_to_script = "../ve/difmap/final_clean_nw"
# Lesha's data
need_downscale_uv = True
template_uvfits_dict = {15.4: "/home/ilya/data/M87Lesha/to_ilya/1228+126.U.2009_05_23C.uvf_cal",
                        8.1: "/home/ilya/data/M87Lesha/to_ilya/1228+126.X.2009_05_23.uvf_cal"}
# MOJAVE data
# need_downscale_uv = False
# template_uvfits_dict = {15.4: "/home/ilya/data/M87Lesha/artMOJAVE/1228+126.u.2020_07_02.uvf",
#                         8.1: "/home/ilya/data/M87Lesha/artMOJAVE/1228+126.x.2006_06_15.uvf"}
##############################################
# No need to change anything below this line #
##############################################
# Plot only jet emission and do not plot counter-jet?
jet_only = True

freq_high = max(freqs_obs_ghz)
freq_low = min(freqs_obs_ghz)

if not only_make_pics:
    # common_beam = create_clean_image_from_fits_file(os.path.join(save_dir, "template_cc_i_{}.fits".format(freq_high))).beam
    for freq in freqs_obs_ghz:
        uvdata = UVData(template_uvfits_dict[freq])
        noise = uvdata.noise(average_freq=False, use_V=False)
        uvdata.zero_data()
        # If one needs to decrease the noise this is the way to do it
        for baseline, baseline_noise_std in noise.items():
            noise.update({baseline: noise_scale_factor*baseline_noise_std})

        jm = JetImage(z=z, n_along=n_along, n_across=n_across,
                      lg_pixel_size_mas_min=lg_pixel_size_mas_min, lg_pixel_size_mas_max=lg_pixel_size_mas_max,
                      jet_side=True, rot=np.deg2rad(-107.0))
        jm.load_image_stokes("I", "{}/{}_jet_image_i_{}.txt".format(jetpol_files_directory, mhdrt_code, 15.4), scale=(15.4/freq)**(-alpha_true))
        # Convert to difmap model format
        jm.save_image_to_difmap_format("{}/true_jet_model_i_{}.txt".format(save_dir, freq))
        # Rotate
        rotate_difmap_model("{}/true_jet_model_i_{}.txt".format(save_dir, freq),
                            "{}/true_jet_model_i_{}_rotated.txt".format(save_dir, freq),
                            PA_deg=107.0)
        # Convolve with beam
        convert_difmap_model_file_to_CCFITS("{}/true_jet_model_i_{}.txt".format(save_dir, freq), "I", common_mapsize,
                                            common_beam, template_uvfits_dict[freq],
                                            "{}/convolved_true_jet_model_i_{}.fits".format(save_dir, freq))
        if jet_only:
            uvdata.substitute([jm])

        else:
            cjm = JetImage(z=z, n_along=n_along, n_across=n_across,
                             lg_pixel_size_mas_min=lg_pixel_size_mas_min, lg_pixel_size_mas_max=lg_pixel_size_mas_max,
                             jet_side=False, rot=np.deg2rad(-107.0))
            cjm.load_image_stokes("I", "{}/{}_cjet_image_i_{}.txt".format(jetpol_files_directory, mhdrt_code, freq))
            jms = TwinJetImage(jm, cjm)
            uvdata.substitute([jms])

        uvdata.noise_add(noise)
        uvdata.save("template_{}.uvf".format(freq), rewrite=True, downscale_by_freq=need_downscale_uv)

        outfname = "model_cc_i_{}.fits".format(freq)
        if os.path.exists(outfname):
            os.unlink(outfname)
        clean_difmap(fname="template_{}.uvf".format(freq),
                     outfname=outfname, outpath=save_dir, stokes="i",
                     mapsize_clean=common_mapsize, path_to_script=path_to_script,
                     show_difmap_output=True,
                     beam_restore=common_beam,
                     dfm_model=os.path.join(save_dir, "model_cc_i_{}.mdl".format(freq)))

# Create image of alpha made from true jet models convolved with beam
itrue_convolved_low = create_image_from_fits_file("{}/convolved_true_jet_model_i_{}.fits".format(save_dir, freq_low))
itrue_convolved_high = create_image_from_fits_file("{}/convolved_true_jet_model_i_{}.fits".format(save_dir, freq_high))
true_convolved_spix_array = np.log(itrue_convolved_low.image/itrue_convolved_high.image)/np.log(freq_low/freq_high)

# Observed images of CLEANed artificial UV-data
ccimages = {freq: create_clean_image_from_fits_file(os.path.join(save_dir, "model_cc_i_{}.fits".format(freq)))
            for freq in freqs_obs_ghz}
ipol = ccimages[freq_low].image
beam = ccimages[freq_low].beam
# Number of pixels in beam
npixels_beam = np.pi*beam[0]*beam[1]/(4*np.log(2)*common_mapsize[1]**2)


std = find_image_std(ipol, beam_npixels=npixels_beam)
print("IPOL image std = {} mJy/beam".format(1000*std))
blc, trc = find_bbox(ipol, level=3*std, min_maxintensity_mjyperbeam=10*std,
                     min_area_pix=10*npixels_beam, delta=10)
if blc[0] == 0: blc = (blc[0]+1, blc[1])
if blc[1] == 0: blc = (blc[0], blc[1]+1)
if trc[0] == ipol.shape: trc = (trc[0]-1, trc[1])
if trc[1] == ipol.shape: trc = (trc[0], trc[1]-1)

# By hand
blc = (459, 474)
trc = (836, 654)

fig = iplot(ccimages[freq_high].image, x=ccimages[freq_high].x, y=ccimages[freq_high].y,
            min_abs_level=3*std, blc=blc, trc=trc, beam=beam, close=True, show_beam=True, show=False,
            contour_color='gray', contour_linewidth=0.25)
fig.savefig(os.path.join(save_dir2, "observed_ipol_{}GHz.png".format(freq_high)), dpi=600, bbox_inches="tight")

fig = iplot(ccimages[freq_low].image, x=ccimages[freq_low].x, y=ccimages[freq_low].y,
            min_abs_level=3*std, blc=blc, trc=trc, beam=beam, close=True, show_beam=True, show=False,
            contour_color='gray', contour_linewidth=0.25)
fig.savefig(os.path.join(save_dir2, "observed_ipol_{}GHz.png".format(freq_low)), dpi=600, bbox_inches="tight")


ipol_arrays = dict()
sigma_ipol_arrays = dict()
masks_dict = dict()
for freq in freqs_obs_ghz:

    ipol = ccimages[freq].image
    ipol_arrays[freq] = ipol

    std = find_image_std(ipol, beam_npixels=npixels_beam)
    masks_dict[freq] = ipol < 3*std
    sigma_ipol_arrays[freq] = np.ones(ipol.shape)*std

common_imask = np.logical_or.reduce([masks_dict[freq] for freq in freqs_obs_ghz])
spix_array = np.log(ipol_arrays[freq_low]/ipol_arrays[freq_high])/np.log(freq_low/freq_high)
# spix_array, sigma_spix_array, spix_chisq_array = spix_map(np.array(freqs_obs_ghz)*1e9, ipol_arrays, sigma_ipol_arrays,
#                                                           mask=common_imask, outfile=None, outdir=None,
#                                                           mask_on_chisq=False, ampcal_uncertainties=None)

fig = iplot(ipol, true_convolved_spix_array, x=ccimages[freq_low].x, y=ccimages[freq_low].y,
            min_abs_level=3*std, colors_mask=common_imask, color_clim=[-0.5, 0.5], blc=blc, trc=trc,
            beam=beam, close=True, colorbar_label=r"$\alpha$", show_beam=True, show=False,
            cmap='bwr', contour_color='black', plot_colorbar=True,
            contour_linewidth=0.25)
fig.savefig(os.path.join(save_dir2, "true_conv_spix_{}GHz_{}GHz.png".format(freq_high, freq_low)), dpi=600, bbox_inches="tight")

fig = iplot(ipol, spix_array, x=ccimages[freq_low].x, y=ccimages[freq_low].y,
            min_abs_level=3*std, colors_mask=common_imask, color_clim=[-1.5, 0.5], blc=blc, trc=trc,
            beam=beam, close=True, colorbar_label=r"$\alpha$", show_beam=True, show=False,
            cmap='gist_ncar', contour_color='black', plot_colorbar=True,
            contour_linewidth=0.25)
fig.savefig(os.path.join(save_dir2, "spix_{}GHz_{}GHz.png".format(freq_high, freq_low)), dpi=600, bbox_inches="tight")