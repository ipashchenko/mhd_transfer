import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
from vlbi_utils import find_image_std, find_bbox, pol_mask, correct_ppol_bias
sys.path.insert(0, '../ve/vlbi_errors')
from image_ops import spix_map, rotm_map, hovatta_find_sigma_pang
from from_fits import create_clean_image_from_fits_file
from image import plot as iplot


# This script generates ROTM and SPIX maps from multifrequency maps obtained from the VLBI data that were generated from
# model brightness distributions.

freqs = [15.4, 12.1, 8.1]
common_mapsize = (512, 0.1)
# Directory with FITS files with CLEAN images made with the same parameters (map, beam)
fits_files_dir = "/home/ilya/data/mf"
save_dir = fits_files_dir
ccimages = dict()
pang_arrays = list()
sigma_pang_arrays = list()
ipol_arrays = list()
sigma_ipol_arrays = list()
beam_npixels = None

stokes = ("i", "q", "u")
for freq in freqs:
    ccimages[freq] = dict()
    for stk in stokes:
        ccimages[freq][stk] = create_clean_image_from_fits_file(os.path.join(fits_files_dir,
                                                                             "model_cc_{}_{}.fits".format(stk, freq)))
        if beam_npixels is None:
            beam = ccimages[freq][stk].beam
            beam_npixels = np.pi*beam[0]*beam[1]/(4*np.log(2)*common_mapsize[1]**2)

    ipol = ccimages[freq]["i"].image
    qpol = ccimages[freq]["q"].image
    upol = ccimages[freq]["u"].image
    ipol_arrays.append(ipol)
    pang_arrays.append(0.5*np.arctan2(upol, qpol))
    sigma_pang_array, sigma_ppol_array = hovatta_find_sigma_pang(ccimages[freq]["q"],
                                                                 ccimages[freq]["u"],
                                                                 ccimages[freq]["i"],
                                                                 sigma_evpa=0, d_term=0, n_ant=1, n_if=1, n_scan=1)
    sigma_pang_arrays.append(sigma_pang_array)
    mask_dict, ppol_quantile = pol_mask({stk.upper(): ccimages[freq][stk].image for stk in stokes}, beam_npixels, n_sigma=3,
                                        return_quantile=True)
    ccimages[freq]["masks"] = mask_dict
    ccimages[freq]["pquantile"] = ppol_quantile

    std = find_image_std(ipol, beam_npixels=beam_npixels)
    sigma_ipol_arrays.append(np.ones(ipol.shape)*std)
    blc, trc = find_bbox(ipol, level=4*std, min_maxintensity_mjyperbeam=10*std,min_area_pix=4*beam_npixels, delta=10)
    if blc[0] == 0: blc = (blc[0]+1, blc[1])
    if blc[1] == 0: blc = (blc[0], blc[1]+1)
    if trc[0] == ipol.shape: trc = (trc[0]-1, trc[1])
    if trc[1] == ipol.shape: trc = (trc[0], trc[1]-1)

    ccimages[freq]["box"] = (blc, trc)
    ccimages[freq]["std"] = std

common_pmask = np.logical_or.reduce([ccimages[freq]["masks"]["P"] for freq in freqs])
common_imask = np.logical_or.reduce([ccimages[freq]["masks"]["I"] for freq in freqs])
rotm_array, sigma_rotm_array, rotm_chisq_array = rotm_map(np.array(freqs)*1e9, pang_arrays, sigma_pang_arrays,
                                                          mask=common_pmask, outfile=None, outdir=None,
                                                          mask_on_chisq=True, plot_pxls=None, outfile_pxls=None)
spix_array, sigma_spix_array, spix_chisq_array = spix_map(np.array(freqs)*1e9, ipol_arrays, sigma_ipol_arrays,
                                                          mask=common_imask, outfile=None, outdir=None,
                                                          mask_on_chisq=False, ampcal_uncertainties=None)

plot_freq = min(freqs)
blc, trc = ccimages[plot_freq]["box"]
fig = iplot(ipol, rotm_array, x=ccimages[plot_freq]["i"].x, y=ccimages[plot_freq]["i"].y,
            min_abs_level=3*ccimages[min(freqs)]["std"], colors_mask=common_pmask, color_clim=None, blc=blc, trc=trc,
            beam=beam, close=True, colorbar_label=r"RM, rad/m$^2$", show_beam=True, show=False,
            cmap='jet', contour_color='black', plot_colorbar=True,
            contour_linewidth=0.25)
fig.savefig(os.path.join(save_dir, "rotm_vlbi.png"), dpi=600, bbox_inches="tight")

fig = iplot(ipol, spix_array, x=ccimages[plot_freq]["i"].x, y=ccimages[plot_freq]["i"].y,
            min_abs_level=3*ccimages[min(freqs)]["std"], colors_mask=common_imask, color_clim=[-1.5, 1], blc=blc, trc=trc,
            beam=beam, close=True, colorbar_label=r"$\alpha$", show_beam=True, show=False,
            cmap='jet', contour_color='black', plot_colorbar=True,
            contour_linewidth=0.25)
fig.savefig(os.path.join(save_dir, "spix_vlbi.png"), dpi=600, bbox_inches="tight")