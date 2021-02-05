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


# Directory to save all
save_dir = "/home/ilya/data/mf"
# Observed frequencies of simulations
freqs_obs_ghz = [8.1, 12.1, 15.4]

# Multiplicative factor for noise added to model visibilities.
noise_scale_factor = 0.1
# Common size of the map and pixel size (mas)
common_mapsize = (1024, 0.1)
# Common beam size (mas, mas, deg)
common_beam = (1.2, 1.2, 0)

jetpol_run_directory = "cmake-build-release"
# C++ code run parameters
z = 0.00436
n_along = 1200
# n_along = 100
n_across = 300
# n_across = 50
lg_pixel_size_mas_min = np.log10(0.02)
# lg_pixel_size_mas_min = np.log10(0.06)
lg_pixel_size_mas_max = np.log10(0.02)
# lg_pixel_size_mas_max = np.log10(0.06)

path_to_script = "../ve/difmap/final_clean_nw"
# Some template UVFITS with full polarization
template_uvfits_dict = {15.4: "1458+718.u.2006_09_06.uvf", 12.1: "1458+718.j.2006_09_06.uvf",
                        8.4: "1458+718.y.2006_09_06.uvf", 8.1: "1458+718.x.2006_09_06.uvf"}

##############################################
# No need to change anything below this line #
##############################################
# Plot only jet emission and do not plot counter-jet?
jet_only = True

for freq in freqs_obs_ghz:
    uvdata = UVData(template_uvfits_dict[freq])
    noise = uvdata.noise(average_freq=False, use_V=False)
    # If one needs to decrease the noise this is the way to do it
    for baseline, baseline_noise_std in noise.items():
        noise.update({baseline: noise_scale_factor*baseline_noise_std})

    stokes = ("I", "Q", "U")
    jms = [JetImage(z=z, n_along=n_along, n_across=n_across,
                    lg_pixel_size_mas_min=lg_pixel_size_mas_min, lg_pixel_size_mas_max=lg_pixel_size_mas_max,
                    jet_side=True) for _ in stokes]
    cjms = [JetImage(z=z, n_along=n_along, n_across=n_across,
                     lg_pixel_size_mas_min=lg_pixel_size_mas_min, lg_pixel_size_mas_max=lg_pixel_size_mas_max,
                     jet_side=False) for _ in stokes]
    for i, stk in enumerate(stokes):
        jms[i].load_image_stokes(stk, "../{}/jet_image_{}_{}.txt".format(jetpol_run_directory, stk.lower(), freq))
        cjms[i].load_image_stokes(stk, "../{}/cjet_image_{}_{}.txt".format(jetpol_run_directory, stk.lower(), freq))

    # List of models (for J & CJ) for all stokes
    js = [TwinJetImage(jms[i], cjms[i]) for i in range(len(stokes))]

    uvdata.zero_data()
    if jet_only:
        uvdata.substitute(jms)
    else:
        uvdata.substitute(js)
    uvdata.noise_add(noise)
    uvdata.save("template_{}.uvf".format(freq), rewrite=True)

    for stk in stokes:
        outfname = "model_cc_{}_{}.fits".format(stk.lower(), freq)
        if os.path.exists(outfname):
            os.unlink(outfname)
        clean_difmap(fname="template_{}.uvf".format(freq),
                     outfname=outfname, outpath=save_dir, stokes=stk.lower(),
                     mapsize_clean=common_mapsize, path_to_script=path_to_script,
                     show_difmap_output=False, beam_restore=common_beam,
                     dfm_model=os.path.join(save_dir, "model_cc_{}_{}.mdl".format(stk.lower(), freq)))
