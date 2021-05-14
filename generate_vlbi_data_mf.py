import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import sys
from jet_image import JetImage, TwinJetImage
sys.path.insert(0, '../ve/vlbi_errors')
from uv_data import UVData
from spydiff import clean_difmap


mhdrt_code = "m1s10g2b123.971372r0.000369_psi_1.000000_dpsi_0.015000"
# Directory to save all
save_dir = "/home/ilya/data/mf/pol"
# Observed frequencies of simulations
# freqs_obs_ghz = [8.1, 12.1, 15.4]
freqs_obs_ghz = [15.4]

# Multiplicative factor for noise added to model visibilities.
noise_scale_factor = 0.1
# Common size of the map and pixel size (mas)
common_mapsize = (1024, 0.1)
# Common beam size (mas, mas, deg)
# common_beam = (1.2, 1.2, 0)
common_beam = None

jetpol_run_directory = save_dir
# C++ code run parameters
n_along = 800
n_across = 150
lg_pixel_size_mas_min = np.log10(0.025)
lg_pixel_size_mas_max = np.log10(0.05)
z = 0.00436

# path_to_script = "/home/ilya/github/stackemall/external_scripts/last/script_clean_rms_last"
path_to_script = "/home/ilya/github/stackemall/external_scripts/last/final_clean_rms_last"
# Some template UVFITS with full polarization
template_uvfits_dict = {15.4: "/home/ilya/data/mf/pol/1228+126.u.2006_06_15.uvf",
                        12.1: "/home/ilya/data/mf/pol/1228+126.j.2006_06_15.uvf",
                        8.4: "/home/ilya/data/mf/pol/1228+126.y.2006_06_15.uvf",
                        8.1: "/home/ilya/data/mf/pol/1228+126.x.2006_06_15.uvf"}
text_boxes = {15.4: "/home/ilya/data/mf/pol/art_wins.txt"}

##############################################
# No need to change anything below this line #
##############################################
# Plot only jet emission and do not plot counter-jet?
jet_only = True

for freq in freqs_obs_ghz:
    print("Converting model image at {} GHz".format(freq))
    uvdata = UVData(template_uvfits_dict[freq])
    noise = uvdata.noise(average_freq=False, use_V=False)
    # If one needs to decrease the noise this is the way to do it
    print("Scaling real noise by coefficient {}".format(noise_scale_factor))
    for baseline, baseline_noise_std in noise.items():
        noise.update({baseline: noise_scale_factor*baseline_noise_std})

    stokes = ("I", "Q", "U")
    # rot=np.deg2rad(-107.0)
    jms = [JetImage(z=z, n_along=n_along, n_across=n_across,
                    lg_pixel_size_mas_min=lg_pixel_size_mas_min, lg_pixel_size_mas_max=lg_pixel_size_mas_max,
                    jet_side=True, rot=np.deg2rad(-107.0)) for _ in stokes]
    # cjms = [JetImage(z=z, n_along=n_along, n_across=n_across,
    #                  lg_pixel_size_mas_min=lg_pixel_size_mas_min, lg_pixel_size_mas_max=lg_pixel_size_mas_max,
    #                  jet_side=False) for _ in stokes]
    for i, stk in enumerate(stokes):
        jms[i].load_image_stokes(stk, "{}/{}_jet_image_{}_{}.txt".format(jetpol_run_directory, mhdrt_code, stk.lower(), freq))
        print("Sum of image intensities: ", np.sum(jms[i].image_intensity()))
        # cjms[i].load_image_stokes(stk, "{}/{}_cjet_image_{}_{}.txt".format(jetpol_run_directory, mhdrt_code, stk.lower(), freq))

    # List of models (for J & CJ) for all stokes
    # js = [TwinJetImage(jms[i], cjms[i]) for i in range(len(stokes))]

    uvdata.zero_data()
    if jet_only:
        print("Substituting data...")
        uvdata.substitute(jms)
    else:
        # uvdata.substitute(js)
        pass
    print("Adding noise...")
    uvdata.rotate_evpa(np.deg2rad(-107.0))
    uvdata.noise_add(noise)
    uvdata.save(os.path.join(save_dir, "template_{}.uvf".format(freq)), rewrite=True)

    for stk in stokes:
        outfname = "model_cc_{}_{}.fits".format(stk.lower(), freq)
        if os.path.exists(outfname):
            os.unlink(outfname)
        print("CLEANing stokes {}...".format(stk))
        clean_difmap(fname="template_{}.uvf".format(freq), path=save_dir,
                     outfname=outfname, outpath=save_dir, stokes=stk.lower(),
                     mapsize_clean=common_mapsize, path_to_script=path_to_script,
                     show_difmap_output=True, beam_restore=common_beam,
                     # text_box=text_boxes[freq],
                     dfm_model=os.path.join(save_dir, "model_cc_{}_{}.mdl".format(stk.lower(), freq)))
