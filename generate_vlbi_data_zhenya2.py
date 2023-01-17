import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import sys
import astropy.io.fits as pf
from vlbi_utils import find_image_std, find_bbox, convert_difmap_model_file_to_CCFITS, CCFITS_to_difmap
sys.path.insert(0, '/home/ilya/github/ve/vlbi_errors')
from uv_data import UVData
from spydiff import clean_difmap
from from_fits import create_clean_image_from_fits_file, create_image_from_fits_file, create_model_from_fits_file
from image_ops import spix_map
from image import plot as iplot
sys.path.insert(0, '/home/ilya/github/alpha')
from alpha_utils import CLEAN_difmap, get_uvrange


# To tweak image parameters
only_make_pics = False

# Directory to save all
save_dir = "/home/ilya/data/M87Lesha/fromZhenya/m87.kq.2018_04_28/art/alpha_paper/uv_clipping"
# Observed frequencies of simulations
freqs_obs_ghz = [24, 43]
freq_high = max(freqs_obs_ghz)
freq_low = min(freqs_obs_ghz)

# Multiplicative factor for noise added to model visibilities.
noise_scale_factor = 1.0
# Common size of the map and pixel size (mas)
common_mapsize = (4096, 0.025)
# Common beam size (mas, mas, deg)
common_beam = (0.42, 0.76, np.deg2rad(-9.8))
common_beam_deg = (0.42, 0.76, -9.8)
# common_beam = (0.5, 0.5, 0)
# common_beam_deg = (0.5, 0.5, 0)

path_to_script = "../ve/difmap/final_clean_nw"
# path_to_script = "/home/ilya/data/M87Lesha/fromZhenya/m87.kq.2018_04_28/art/script_clean_rms_last"
template_uvfits_dict = {24: "/home/ilya/data/M87Lesha/fromZhenya/m87.kq.2018_04_28/m87.k.2018_04_28.fin_uvf_cal",
                        43: "/home/ilya/data/M87Lesha/fromZhenya/m87.kq.2018_04_28/m87.q.2018_04_28.fin_uvf_cal"}
uvrange_low = get_uvrange(template_uvfits_dict[24])
uvrange_high = get_uvrange(template_uvfits_dict[43])
uvrange = (uvrange_high[0], uvrange_low[1])
print("Common uvrage = ", uvrange)
# model_ccfits_dict = {24: "/home/ilya/data/M87Lesha/fromZhenya/m87.kq.2018_04_28/art/alpha_paper/uv_clipping/spi.ak.circ05.px002.ifits",
#                      43: "/home/ilya/data/M87Lesha/fromZhenya/m87.kq.2018_04_28/art/alpha_paper/uv_clipping/spi.aq.circ05.px002.ifits"}
model_ccfits_dict = {24: "/home/ilya/data/M87Lesha/fromZhenya/m87.kq.2018_04_28/m87.k.2018_04_28.fin_ifits",
                     43: "/home/ilya/data/M87Lesha/fromZhenya/m87.kq.2018_04_28/m87.q.2018_04_28.fin_ifits"}
text_boxes = {24: "/home/ilya/data/M87Lesha/fromZhenya/m87.kq.2018_04_28/m87.k.2018_04_28.fin_win",
              43: "/home/ilya/data/M87Lesha/fromZhenya/m87.kq.2018_04_28/m87.q.2018_04_28.fin_wins"}

original_cc = create_clean_image_from_fits_file(model_ccfits_dict[24])
# common_beam = original_cc.beam
# common_beam = (common_beam[0], common_beam[1], np.rad2deg(common_beam[2]))
print("===================")
print("Common beam will be : ", common_beam)
print("===================")

# Create ground truths
if not only_make_pics:
    for freq in freqs_obs_ghz:
        CCFITS_to_difmap(model_ccfits_dict[freq], os.path.join(save_dir, "cc_true_{}.mdl".format(freq)))
        convert_difmap_model_file_to_CCFITS(os.path.join(save_dir, "cc_true_{}.mdl".format(freq)), "I",
                                            common_mapsize, common_beam_deg, template_uvfits_dict[freq],
                                            os.path.join(save_dir, "ccmap_true_{}.fits".format(freq)))
itrue_high = pf.getdata(os.path.join(save_dir, "ccmap_true_{}.fits".format(freq_high))).squeeze()
itrue_low = pf.getdata(os.path.join(save_dir, "ccmap_true_{}.fits".format(freq_low))).squeeze()
true_spix_array = np.log(itrue_low/itrue_high)/np.log(freq_low/freq_high)

# if not only_make_pics:
#     for freq in freqs_obs_ghz:
#         # TODO: Use new clean function with boxes!
#         clean_difmap(fname=template_uvfits_dict[freq], path=save_dir,
#                      outfname="cc_common_original_{}.fits".format(freq), outpath=save_dir, stokes="i",
#                      mapsize_clean=common_mapsize, path_to_script=path_to_script,
#                      show_difmap_output=True,
#                      beam_restore=common_beam, text_box=text_boxes[freq],
#                      dfm_model=os.path.join(save_dir, "model_cc_common_original_i_{}.mdl".format(freq)))
#
#
# itrue_high_image = create_image_from_fits_file(os.path.join(save_dir, "cc_common_original_{}.fits".format(freq_high))).image
# itrue_low_image = create_image_from_fits_file(os.path.join(save_dir, "cc_common_original_{}.fits".format(freq_low))).image
# true_spix_array_image = np.log(itrue_low_image/itrue_high_image)/np.log(freq_low/freq_high)

if not only_make_pics:
    # common_beam = create_clean_image_from_fits_file(os.path.join(save_dir, "template_cc_i_{}.fits".format(freq_high))).beam
    for freq in freqs_obs_ghz:
        uvdata = UVData(template_uvfits_dict[freq])
        noise = uvdata.noise(average_freq=False, use_V=False)
        uvdata.zero_data()
        # If one needs to decrease the noise this is the way to do it
        for baseline, baseline_noise_std in noise.items():
            noise.update({baseline: noise_scale_factor*baseline_noise_std})

        jm = create_model_from_fits_file(model_ccfits_dict[freq])
        uvdata.substitute([jm])

        uvdata.noise_add(noise)
        uvdata.save(os.path.join(save_dir, "template_{}.uvf".format(freq)), rewrite=True, downscale_by_freq=True)

        outfname = "model_cc_i_{}.fits".format(freq)
        if os.path.exists(outfname):
            os.unlink(outfname)
        # clean_difmap(fname="template_{}.uvf".format(freq), path=save_dir,
        #              outfname=outfname, outpath=save_dir, stokes="i",
        #              mapsize_clean=common_mapsize, path_to_script=path_to_script,
        #              show_difmap_output=True,
        #              beam_restore=common_beam,
        #              text_box=text_boxes[freq],
        #              dfm_model=os.path.join(save_dir, "model_cc_i_{}.mdl".format(freq)))

        CLEAN_difmap("template_{}.uvf".format(freq), "i", common_mapsize, "model_cc_i_{}".format(freq), restore_beam=common_beam_deg,
                     boxfile=text_boxes[freq], working_dir=save_dir,
                     uvrange=uvrange,
                     box_clean_nw_niter=1000, clean_gain=0.03, dynam_su=20, dynam_u=6, deep_factor=1.0,
                     remove_difmap_logs=False, save_noresid=False, save_resid_only=False, save_dfm=False,
                     noise_to_use="F")


ccimages = {freq: create_image_from_fits_file(os.path.join(save_dir, "model_cc_i_{}.fits".format(freq)))
            for freq in freqs_obs_ghz}
ipol = ccimages[freq_low].image
beam = common_beam
# Number of pixels in beam
npixels_beam = np.pi*beam[0]*beam[1]/(4*np.log(2)*common_mapsize[1]**2)


std = find_image_std(ipol, beam_npixels=npixels_beam)
print("IPOL image std = {} mJy/beam".format(1000*std))
# blc, trc = find_bbox(ipol, level=3*std, min_maxintensity_mjyperbeam=10*std,
#                      min_area_pix=10*npixels_beam, delta=10)
# if blc[0] == 0: blc = (blc[0]+1, blc[1])
# if blc[1] == 0: blc = (blc[0], blc[1]+1)
# if trc[0] == ipol.shape: trc = (trc[0]-1, trc[1])
# if trc[1] == ipol.shape: trc = (trc[0], trc[1]-1)


# By hand
# blc = (1950, 1950)
blc = (1960, 1960)
# trc = (2250, 2150)
trc = (2300, 2180)

# IPOL contours at high freq
fig = iplot(ipol, x=ccimages[freq_high].x, y=ccimages[freq_high].y,
            min_abs_level=3*std, blc=blc, trc=trc, beam=beam, close=True, show_beam=True, show=False,
            contour_color='gray', contour_linewidth=0.25)
fig.savefig(os.path.join(save_dir, "observed_ipol_{}GHz.png".format(freq_high)), dpi=600, bbox_inches="tight")

# ALPHA
# fig = iplot(ipol, fpol, x=ccimages["I"].x, y=ccimages["I"].y,
#             min_abs_level=3*std, colors_mask=masks_dict["P"], color_clim=[0, 0.7], blc=blc, trc=trc,
#             beam=beam, close=True, colorbar_label="m", show_beam=True, show=False,
#             cmap='gnuplot', contour_color='black', plot_colorbar=True,
#             contour_linewidth=0.25)
# fig.savefig(os.path.join(save_dir, "observed_fpol.png"), dpi=600, bbox_inches="tight")


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
bias_spix_array = spix_array - true_spix_array
# spix_array, sigma_spix_array, spix_chisq_array = spix_map(np.array(freqs_obs_ghz)*1e9, ipol_arrays, sigma_ipol_arrays,
#                                                           mask=common_imask, outfile=None, outdir=None,
#                                                           mask_on_chisq=False, ampcal_uncertainties=None)

fig = iplot(ipol, spix_array, x=ccimages[freq_low].x, y=ccimages[freq_low].y,
            min_abs_level=5*std, colors_mask=common_imask, color_clim=[-1.5, 0.5], blc=blc, trc=trc,
            beam=beam, close=True, colorbar_label=r"$\alpha$", show_beam=True, show=False,
            cmap='nipy_spectral', contour_color='black', plot_colorbar=True,
            contour_linewidth=0.25)
fig.savefig(os.path.join(save_dir, "spix_observed_{}GHz_{}GHz.png".format(freq_high, freq_low)), dpi=600, bbox_inches="tight")

fig = iplot(ipol, bias_spix_array, x=ccimages[freq_low].x, y=ccimages[freq_low].y,
            min_abs_level=5*std, colors_mask=common_imask, color_clim=[-0.21, 0.21], blc=blc, trc=trc,
            beam=beam, close=True, colorbar_label=r"$\alpha$ bias", show_beam=True, show=False,
            cmap='coolwarm', contour_color='black', plot_colorbar=True,
            contour_linewidth=0.25, show_xlabel_on_current_axes=True, show_ylabel_on_current_axes=True)
fig.savefig(os.path.join(save_dir, "spix_bias_{}GHz_{}GHz.png".format(freq_high, freq_low)), dpi=600, bbox_inches="tight")

# fig = iplot(ipol, spix_array - bias_spix_array, x=ccimages[freq_low].x, y=ccimages[freq_low].y,
#             min_abs_level=5*std, colors_mask=common_imask, color_clim=[-1.5, 0.5], blc=blc, trc=trc,
#             beam=beam, close=True, colorbar_label=r"$\alpha$", show_beam=True, show=False,
#             cmap='nipy_spectral', contour_color='black', plot_colorbar=True,
#             contour_linewidth=0.25)
# fig.savefig(os.path.join(save_dir, "spix_bias_corrected_{}GHz_{}GHz.png".format(freq_high, freq_low)), dpi=600, bbox_inches="tight")
#
# fig = iplot(ipol, spix_array - 2*bias_spix_array, x=ccimages[freq_low].x, y=ccimages[freq_low].y,
#             min_abs_level=5*std, colors_mask=common_imask, color_clim=[-1.5, 0.5], blc=blc, trc=trc,
#             beam=beam, close=True, colorbar_label=r"$\alpha$", show_beam=True, show=False,
#             cmap='nipy_spectral', contour_color='black', plot_colorbar=True,
#             contour_linewidth=0.25)
# fig.savefig(os.path.join(save_dir, "spix_bias_corrected_x2_{}GHz_{}GHz.png".format(freq_high, freq_low)), dpi=600, bbox_inches="tight")
#
# fig = iplot(ipol, true_spix_array_image - bias_spix_array, x=ccimages[freq_low].x, y=ccimages[freq_low].y,
#             min_abs_level=5*std, colors_mask=common_imask, color_clim=[-1.5, 0.5], blc=blc, trc=trc,
#             beam=beam, close=True, colorbar_label=r"$\alpha$", show_beam=True, show=False,
#             cmap='nipy_spectral', contour_color='black', plot_colorbar=True,
#             contour_linewidth=0.25)
# fig.savefig(os.path.join(save_dir, "true_spix_image_bias_corrected_{}GHz_{}GHz.png".format(freq_high, freq_low)), dpi=600, bbox_inches="tight")
#
# fig = iplot(ipol, true_spix_array, x=ccimages[freq_low].x, y=ccimages[freq_low].y,
#             min_abs_level=5*std, colors_mask=common_imask, color_clim=[-1.5, 0.5], blc=blc, trc=trc,
#             beam=beam, close=True, colorbar_label=r"$\alpha$", show_beam=True, show=False,
#             cmap='nipy_spectral', contour_color='black', plot_colorbar=True,
#             contour_linewidth=0.25)
# fig.savefig(os.path.join(save_dir, "spix_true_{}GHz_{}GHz.png".format(freq_high, freq_low)), dpi=600, bbox_inches="tight")