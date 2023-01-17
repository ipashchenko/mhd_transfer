from abc import ABC
import os
import functools
from astropy import units as u, cosmology
import numpy as np
import datetime
import matplotlib
label_size = 16
matplotlib.rcParams['ytick.labelsize'] = label_size
matplotlib.rcParams['axes.titlesize'] = label_size
matplotlib.rcParams['axes.labelsize'] = label_size
matplotlib.rcParams['font.size'] = label_size
matplotlib.rcParams['legend.fontsize'] = label_size
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, LogLocator
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
try:
    from fourier import FINUFFT_NUNU
except ImportError:
    raise Exception("Install pynufft")
    FINUFFT_NUNU = None



def convert_difmap_model_file_to_CCFITS(difmap_model_file, stokes, mapsize, restore_beam, uvfits_template, out_ccfits,
                                        show_difmap_output=True):
    """
    Using difmap-formated model file (e.g. flux, r, theta) obtain convolution of your model with the specified beam.

    :param difmap_model_file:
        Difmap-formated model file. Use ``JetImage.save_image_to_difmap_format`` to obtain it.
    :param stokes:
        Stokes parameter.
    :param mapsize:
        Iterable of image size and pixel size (mas).
    :param restore_beam:
        Beam to restore: bmaj(mas), bmin(mas), bpa(deg).
    :param uvfits_template:
        Template uvfits observation to use. Difmap can't read model without having observation at hand.
    :param out_ccfits:
        File name to save resulting convolved map.
    :param show_difmap_output: (optional)
        Boolean. Show Difmap output? (default: ``True``)
    """
    stamp = datetime.datetime.now()
    command_file = "difmap_commands_{}".format(stamp.isoformat())
    difmapout = open(command_file, "w")
    difmapout.write("observe " + uvfits_template + "\n")
    difmapout.write("select " + stokes + "\n")
    difmapout.write("rmodel " + difmap_model_file + "\n")
    difmapout.write("mapsize " + str(mapsize[0]) + "," + str(mapsize[1]) + "\n")
    print("Restoring difmap model with BEAM : bmin = " + str(restore_beam[1]) + ", bmaj = " + str(restore_beam[0]) + ", " + str(restore_beam[2]) + " deg")
    # default dimfap: false,true (parameters: omit_residuals, do_smooth)
    difmapout.write("restore " + str(restore_beam[1]) + "," + str(restore_beam[0]) + "," + str(restore_beam[2]) +
                    "," + "true,false" + "\n")
    difmapout.write("wmap " + out_ccfits + "\n")
    difmapout.write("exit\n")
    difmapout.close()

    shell_command = "difmap < " + command_file + " 2>&1"
    if not show_difmap_output:
        shell_command += " >/dev/null"
    os.system(shell_command)
    os.unlink(command_file)


class JetImage(ABC):
    cosmo = cosmology.WMAP9

    def __init__(self, z, n_along, n_across, lg_pixel_size_mas_min, lg_pixel_size_mas_max, ft_class=FINUFFT_NUNU,
                 jet_side=True, dx=0.0, dy=0.0, rot=0.0):
        self.jet_side = jet_side
        self.z = z
        self.lg_pixel_size_mas_min = lg_pixel_size_mas_min
        self.lg_pixel_size_mas_max = lg_pixel_size_mas_max
        self.n_along = n_along
        self.n_across = n_across
        resolutions = np.logspace(lg_pixel_size_mas_min, lg_pixel_size_mas_max, n_along)
        # 2D array of u.angle pixsizes
        self.pixsize_array = np.tile(resolutions, n_across).reshape(n_across, n_along).T*u.mas
        self.ft_class = ft_class
        self.ft_instance = None
        self.calculate_grid()
        # self.stokes_dict = {stokes: None for stokes in ("I", "Q", "U", "V")}
        self._image = None
        self._image_tau = None
        self._image_alpha = None
        self.stokes = None

        # For compatibility with FT realization and possible future use
        # shift in mas
        self.dx = dx
        self.dy = dy
        # rotation angle in rad
        self.rot = rot

    def plot_resolutions(self):
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 1)
        axes.plot(np.cumsum(self.pixsize_array[:, 0]),
                  np.cumsum(self.pixsize_array[:, :int(self.n_across/2)], axis=1)[:, -1])
        axes.set_xlabel("Along distance, mas")
        axes.set_ylabel("Across distance, mas")
        axes.set_aspect("auto")
        return fig

    def halfphi_app_max(self):
        """
        Maximum half opening-angle [rad] that can be imaged with given resolution.
        """
        return np.arctan(np.sum(self.pixsize_array[-1, :int(self.n_across/2)]) / np.sum(self.pixsize_array[:, 0]))

    def calculate_grid(self):
        """
        Calculate grid of ``(r_ob, d)`` - each point is in the center of the
        corresponding pixel. Thus, ``JetModelZoom.ft`` method should not shift
        phases on a half of a pixel. ``r_ob`` and ``d`` are in parsecs.
        """
        pc_x = np.cumsum(self.pix_to_pc, axis=0)-self.pix_to_pc/2
        pc_y_up = self.pix_to_pc[:, :self.pix_to_pc.shape[1]//2][::-1]
        pc_y_low = self.pix_to_pc[:, self.pix_to_pc.shape[1]//2:]
        pc_y_up = (np.cumsum(pc_y_up, axis=1) - pc_y_up/2)[::-1]
        pc_y_low = np.cumsum(pc_y_low, axis=1) - pc_y_low/2
        pc_y = np.hstack((pc_y_up[:, ::-1], -pc_y_low))
        self.r_ob = pc_x
        if not self.jet_side:
            self.r_ob *= -1
        # FIXME: In analogy with older convention we need ``-`` here
        self.d = -pc_y

    @property
    def pc_to_mas(self):
        return (u.pc/self.ang_to_dist).to(u.mas)

    @property
    def r_ob_mas(self):
        return self.r_ob*self.pc_to_mas

    @property
    def d_mas(self):
        return self.d*self.pc_to_mas

    @property
    def imgsize(self):
        return np.max(np.cumsum(self.pixsize_array, axis=0)), \
               np.max(np.cumsum(self.pixsize_array, axis=1))

    @property
    def img_extent(self):
        s0, s1 = self.imgsize[0].value, self.imgsize[1].value
        if self.jet_side:
            return 0, s0, -s1/2, s1/2
        else:
            return -s0, 0, -s1/2, s1/2

    def load_image_stokes(self, stokes, image_stokes_file, scale=1.0):
        self.stokes = stokes
        image = np.loadtxt(image_stokes_file)
        print("Loaded image with total flux = {} Jy, shape = {}".format(np.nansum(image), image.shape))
        print("Scaling image x ", scale)
        image *= scale
        image[np.isnan(image)] = 0.0
        self._image = image

    def load_image_tau(self, image_tau_file):
        self._image_tau = np.loadtxt(image_tau_file)

    def load_image_alpha(self, image_alpha):
        self._image_alpha = image_alpha

    def save_image_to_difmap_format(self, difmap_format_file, scale=1.0):
        with open(difmap_format_file, "w") as fo:
            for idx, imval in np.ndenumerate(self._image):
                if imval == 0:
                    continue
                j, i = idx
                dec = -self.r_ob_mas[i, j].value
                ra = -self.d_mas[i, j].value
                # print("RA = {}, DEC = {}".format(ra, dec))
                fo.write("{} {} {}\n".format(imval*scale, np.hypot(ra, dec), np.rad2deg(np.arctan2(ra, dec))))

    def image(self):
        return self._image
    
    def image_intensity(self):
        # Factor that accounts non-uniform pixel size in plotting
        factor = (self.pixsize_array/np.min(self.pixsize_array))**2
        return self.image()/factor.T

    def image_tau(self):
        return self._image_tau

    def ft(self, uv):
        mas_to_rad = u.mas.to(u.rad)
        rot = np.array([[np.cos(self.rot), -np.sin(self.rot)],
                        [np.sin(self.rot), np.cos(self.rot)]])

        # No need in half pixel shifts cause coordinates are already at pixel
        # centers
        shift = [self.dx*mas_to_rad, self.dy*mas_to_rad]
        result = np.exp(-2.0*np.pi*1j*(uv @ shift))
        uv = uv @ rot

        x = (self.d*u.pc/self.ang_to_dist).to(u.rad).value
        y = (self.r_ob*u.pc/self.ang_to_dist).to(u.rad).value
        ft_instance = self.ft_class(uv, x.ravel(), y.ravel())
        img = self.image()
        result *= ft_instance.forward(img.T.ravel())
        del ft_instance, x, y

        return result

    @property
    @functools.lru_cache()
    def ang_to_dist(self):
        return self.cosmo.kpc_proper_per_arcmin(self.z)

    @property
    @functools.lru_cache()
    def pix_to_pc(self):
        """
        2D array of pixel sizes in parsecs.
        """
        return (self.pixsize_array * self.ang_to_dist).to(u.pc).value

    def plot_contours(self, nlevels=25, zoom_fr=1.0, outfile=None, frac_min=0.0001, axis_units="mas",
                      count_levels_from_image_min=True, loglevs=True, logstep=np.sqrt(2), contour_cmap="viridis"):
        assert axis_units in ("pc", "mas")
        # zoom fraction - fraction of the original image to show. 0.5 means that show only first half of the image
        assert 0.0 < zoom_fr <= 1.0
        cumsum_length = np.cumsum(self.pixsize_array[:, 0]).value
        original_length = cumsum_length[-1]
        length_to_show = zoom_fr*original_length
        index_to_show = np.argmin(np.abs(cumsum_length - length_to_show)) + 1

        image = self.image_intensity()
        fig, axes = plt.subplots(1, 1)
        image = image[:, :index_to_show]
        image = 1e6*np.ma.array(image, mask=image == 0)
        if count_levels_from_image_min:
            count_from = image.min()
        else:
            count_from = 0.0
        if loglevs:
            levels = LogLocator(base=logstep).tick_values(count_from+frac_min*image.max(), image.max())
            norm = LogNorm(vmin=image.min(), vmax=image.max())
        else:
            levels = MaxNLocator(nbins=nlevels).tick_values(count_from+frac_min*image.max(), image.max())
            norm = None
        print(levels)
        # Contours are *point* based plots, so it is suitable for ``d`` and
        # ``r_ob`` that are centers of pixels.
        if axis_units == "mas":
            x = self.r_ob_mas[:index_to_show, :]
            y =self.d_mas[:index_to_show, :]
            axes.set_ylabel("DEC, mas")
            axes.set_xlabel("RA, mas")
        else:
            x = self.r_ob[:index_to_show, :]
            y =self.d[:index_to_show, :]
            axes.set_ylabel(r"$d$, pc")
            axes.set_xlabel(r"$r_{\rm ob}$, pc")

        cf = axes.contour(x, y, image.T, levels=levels, cmap=contour_cmap, alpha=0.5, norm=norm)
        axes.set_aspect("equal")

        # Make a colorbar with label and units
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("right", size="10%", pad=0.00)
        cb = fig.colorbar(cf, cax=cax)
        # Intensity is in Jy per minimal pixel
        cb.set_label(r"$Intensity, \mu Jy/pix$")
        if outfile:
            fig.savefig(outfile, dpi=300, bbox_inches="tight")
        return fig

    def plot(self, save_picture_file_name=None, save_dir=None, aspect="equal", Nan2zero=True, log=True, axis_units="mas", zoom_fr=1.0, cmap="magma",
             figsize=None, cblabel="Intensity, Jy/pix", scale_by_pixsize=True, beta_app_max=None, beta_app_min=None,
             plot_polarization=False, model_txt_files_dict=None):

        if save_dir is None:
            save_dir = os.getcwd()

        # Factor that accounts non-uniform pixel size in plotting
        assert axis_units in ("pc", "mas")
        factor = (self.pixsize_array/np.min(self.pixsize_array))**2
        factor = factor.value
        if scale_by_pixsize:
            image = self.image()/factor.T
        else:
            image = self.image()
        # image[image < 1e-14] = 0
        min_positive = np.min(image[image > 0])
        print("Min positive = ", min_positive)
        from scipy.stats import scoreatpercentile
        first_contour = scoreatpercentile(image[image > min_positive], 0.1)
        # print(first_contour)
        # if log and Nan2zero:
        #     image[image == 0.0] = 1e-12

        # zoom fraction - fraction of the original image to show. 0.5 means that show only first half of the image
        assert 0.0 < zoom_fr <= 1.0
        cumsum_length = np.cumsum(self.pixsize_array[:, 0]).value
        original_length = cumsum_length[-1]
        length_to_show = zoom_fr*original_length
        index_to_show = np.argmin(np.abs(cumsum_length - length_to_show)) + 1

        fig, axes = plt.subplots(1, 1, figsize=figsize)
        cmap = plt.get_cmap(cmap)
        if log:
            # norm = LogNorm(vmin=min_positive, vmax=image.max())
            norm = LogNorm(vmin=first_contour, vmax=0.33*image.max())
        else:
            norm = None

        # Here X and Y are 2D arrays of bounds, so ``image`` should be the value
        # *inside* those bounds. Therefore, we should remove the last value from
        # the ``image`` array. Currently we are not doing it.
        image = image[:, :index_to_show]
        if axis_units == "mas":
            x = self.r_ob_mas[:index_to_show, :].value
            y = self.d_mas[:index_to_show, :].value
            axes.set_ylabel("DEC, mas")
            axes.set_xlabel("RA, mas")
        else:
            x = self.r_ob[:index_to_show, :]
            y = self.d[:index_to_show, :]
            axes.set_ylabel(r"$d$, pc")
            axes.set_xlabel(r"$r_{\rm ob}$, pc")

        im = axes.pcolormesh(x, y, image.T, norm=norm, cmap=cmap, shading="gouraud", vmin=beta_app_min, vmax=beta_app_max)
        axes.set_aspect(aspect)

        # Make a colorbar with label and units
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("right", size="10%", pad=0.00)
        cb = fig.colorbar(im, cax=cax)
        # Intensity is in Jy per minimal pixel
        cb.set_label(cblabel)

        if save_picture_file_name:
            fig.savefig(os.path.join(save_dir, "I_" + save_picture_file_name), dpi=600, bbox_inches="tight")


        if plot_polarization:
            image_i = image.T
            if model_txt_files_dict is None:
                raise Exception("Provide model_txt_files_dict for loading Stokes Q, U, V!")
            jm.load_image_stokes("Q", model_txt_files_dict["Q"])
            if scale_by_pixsize:
                image = self.image()/factor.T
            else:
                image = self.image()
            image_q = image

            jm.load_image_stokes("U", model_txt_files_dict["U"])
            if scale_by_pixsize:
                image = self.image()/factor.T
            else:
                image = self.image()
            image_u = image

            image_p = np.hypot(image_q, image_u)
            image_chi = 0.5*np.arctan2(image_u, image_q)


            fig, axes = plt.subplots(1, 1, figsize=figsize)
            cmap = plt.get_cmap(cmap)
            if log:
                # norm = LogNorm(vmin=min_positive, vmax=image.max())
                norm = LogNorm(vmin=first_contour, vmax=0.33*image.max())
            else:
                norm = None

            # Here X and Y are 2D arrays of bounds, so ``image`` should be the value
            # *inside* those bounds. Therefore, we should remove the last value from
            # the ``image`` array. Currently we are not doing it.
            image = image[:, :index_to_show]
            if axis_units == "mas":
                x = self.r_ob_mas[:index_to_show, :].value
                y = self.d_mas[:index_to_show, :].value
                axes.set_ylabel("DEC, mas")
                axes.set_xlabel("RA, mas")
            else:
                x = self.r_ob[:index_to_show, :]
                y = self.d[:index_to_show, :]
                axes.set_ylabel(r"$d$, pc")
                axes.set_xlabel(r"$r_{\rm ob}$, pc")

            im = axes.pcolormesh(x, y, image_p.T, norm=norm, cmap=cmap, shading="gouraud", vmin=beta_app_min, vmax=beta_app_max)
            axes.set_aspect(aspect)

            # Make a colorbar with label and units
            divider = make_axes_locatable(axes)
            cax = divider.append_axes("right", size="10%", pad=0.00)
            cb = fig.colorbar(im, cax=cax)
            # Intensity is in Jy per minimal pixel
            cb.set_label("Pol. Intensity, Jy/pix")

            fig.savefig(os.path.join(save_dir, "P_" + save_picture_file_name), dpi=600, bbox_inches="tight")

        return fig

    def plot_alpha(self, nlevels=25, zoom_fr=1.0, outfile=None, frac_min=0.0001, axis_units="mas", figsize=None,
                   alpha_min=None, alpha_max=None,
                   count_levels_from_image_min=True, loglevs=True, logstep=np.sqrt(2)):
        assert axis_units in ("pc", "mas")
        # zoom fraction - fraction of the original image to show. 0.5 means that show only first half of the image
        assert 0.0 < zoom_fr <= 1.0
        cumsum_length = np.cumsum(self.pixsize_array[:, 0]).value
        original_length = cumsum_length[-1]
        length_to_show = zoom_fr*original_length
        index_to_show = np.argmin(np.abs(cumsum_length - length_to_show)) + 1

        image = self._image_alpha
        if image is None:
            raise Exception("Load alpha image first!")
        fig, axes = plt.subplots(1, 1, figsize=figsize)
        image = image[:, :index_to_show]
        if count_levels_from_image_min:
            count_from = image.min()
        else:
            count_from = alpha_min
        levels = MaxNLocator(nbins=nlevels).tick_values(count_from, image.max())
        norm = None
        # Contours are *point* based plots, so it is suitable for ``d`` and
        # ``r_ob`` that are centers of pixels.
        if axis_units == "mas":
            x = self.r_ob_mas[:index_to_show, :].value
            y =self.d_mas[:index_to_show, :].value
            axes.set_ylabel("DEC, mas")
            axes.set_xlabel("RA, mas")
        else:
            x = self.r_ob[:index_to_show, :]
            y =self.d[:index_to_show, :]
            axes.set_ylabel(r"$d$, pc")
            axes.set_xlabel(r"$r_{\rm ob}$, pc")

        cf = axes.contour(x, y, image.T, levels=levels, colors="gray", alpha=0.5, norm=norm)
        im = axes.pcolormesh(x, y, image[:, :index_to_show].T, vmin=alpha_min, vmax=alpha_max, cmap="jet", shading="gouraud")
        axes.set_aspect("equal")

        # Make a colorbar with label and units
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("right", size="10%", pad=0.00)
        cb = fig.colorbar(im, cax=cax)
        # Intensity is in Jy per minimal pixel
        cb.set_label(r"$\alpha$")
        if outfile:
            fig.savefig(outfile, dpi=300, bbox_inches="tight")
        return fig


class TwinJetImage(object):
    """
    Wrapper class used for FT of jet and counter-jet.
    """
    def __init__(self, jet, counterjet):
        self.jet = jet
        self.counterjet = counterjet
        assert jet.stokes == counterjet.stokes
        self.stokes = jet.stokes

    def ft(self, uv):
        return self.jet.ft(uv) + self.counterjet.ft(uv)

    def plot_contours(self, nlevels=25, zoom_fr=1.0, outfile=None, aspect="equal"):
        # zoom fraction - fraction of the original image to show. 0.5 means that show only first half of the image
        assert 0.0 < zoom_fr <= 1.0
        cumsum_length = np.cumsum(self.jet.pixsize_array[:, 0]).value
        original_length = cumsum_length[-1]
        length_to_show = zoom_fr*original_length
        index_to_show = np.argmin(np.abs(cumsum_length - length_to_show)) + 1

        # (2000, 200)
        image = np.hstack((self.counterjet.image_intensity()[:, :index_to_show][::-1],
                           self.jet.image_intensity()[:, :index_to_show]))
        image = 1e6*np.ma.array(image, mask=image == 0)
        r_ob = np.vstack((self.counterjet.r_ob[:index_to_show, :],
                          self.jet.r_ob[:index_to_show, :]))
        d = np.vstack((self.counterjet.d[:index_to_show, :],
                       self.jet.d[:index_to_show, :]))

        import matplotlib.pyplot as plt
        from matplotlib.ticker import MaxNLocator
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        fig, axes = plt.subplots(1, 1)
        levels = MaxNLocator(nbins=nlevels).tick_values(image.min()+0.001*image.max(),
                                                        image.max())
        # Contours are *point* based plots, so it is suitable for ``d`` and
        # ``r_ob`` that are centers of pixels.
        cf = axes.contour(r_ob, d, image.T, levels=levels, colors="black")
        axes.set_ylabel(r"$d$, pc")
        axes.set_xlabel(r"$r_{\rm ob}$, pc")
        axes.set_aspect(aspect)

        # Make a colorbar with label and units
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("right", size="5%", pad=0.00)
        cb = fig.colorbar(cf, cax=cax)
        # Intensity is in Jy per minimal pixel
        cb.set_label(r"$Intensity, \mu Jy/pix$")
        if outfile:
            fig.savefig(outfile, dpi=300, bbox_inches="tight")
        return fig


def plot_interpolated(mhd_code):

    Psi = np.loadtxt("{}_Psi_interpolated.txt".format(mhd_code))
    Psi_m = np.ma.array(Psi, mask = Psi == 1)
    plt.matshow(Psi_m[::,:125].T[::-1, ::], cmap="hsv", aspect="auto");plt.colorbar();plt.show()

    jsq = np.loadtxt("{}_jsq_interpolated.txt".format(mhd_code))
    jsq_m = np.ma.array(jsq, mask=jsq == 0)
    plt.matshow(np.log(jsq_m[::, :125].T[::-1, ::]), cmap="hsv", aspect="auto");plt.colorbar();plt.show()

    bphi = np.loadtxt("{}_beta_phi_interpolated.txt".format(mhd_code))
    bphi_m = np.ma.array(bphi, mask = bphi == 0)
    plt.matshow(bphi_m[::,:125].T[::-1, ::], cmap="hsv", aspect="auto");plt.colorbar();plt.show()

    Gamma = np.loadtxt("{}_Gamma_interpolated.txt".format(mhd_code))
    Gamma_m = np.ma.array(Gamma, mask = Gamma == 1)
    plt.matshow(Gamma_m[::,:125].T[::-1, ::], cmap="hsv", aspect="auto");plt.colorbar();plt.show()

    N = np.loadtxt("{}_N_interpolated.txt".format(mhd_code))
    N_m = np.ma.array(N, mask = N == 0)
    plt.matshow(np.log(N_m[::,:125].T[::-1, ::]), cmap="hsv", aspect="auto");plt.colorbar();plt.show()

    B_phi = np.loadtxt("{}_B_phi_interpolated.txt".format(mhd_code))
    B_phi_m = np.ma.array(B_phi, mask = B_phi == 0)
    plt.matshow(np.log(B_phi_m[::,:125].T[::-1, ::]), cmap="hsv", aspect="auto");plt.colorbar();plt.show()

    B_p = np.loadtxt("{}_B_p_interpolated.txt".format(mhd_code))
    B_p_m = np.ma.array(B_p, mask = B_p == 0)
    plt.matshow(np.log(B_p_m[::,:125].T[::-1, ::]), cmap="hsv", aspect="auto");plt.colorbar();plt.show()

    B_tot = np.hypot(B_p_m, B_phi_m/Gamma_m)
    plt.matshow(np.log((B_tot**2)[::,:125].T[::-1, ::]), cmap="hsv", aspect="auto");plt.colorbar();plt.show()

    beta = np.sqrt(Gamma_m**2-1)/Gamma_m
    D = 1/(Gamma_m*(1-beta*np.cos(np.deg2rad(17))))
    boost = D**2.5
    plt.matshow(boost[::,:125].T[::-1, ::], cmap="hsv", aspect="auto");plt.colorbar();plt.show()

    # n_plasma
    plt.matshow(np.log10(N_m*boost/Gamma_m)[::,:125].T[::-1, ::], cmap="hsv", aspect="auto");plt.colorbar();plt.show()
    # B_tot
    plt.matshow(np.log10(B_tot*boost)[::,:125].T[::-1, ::], cmap="hsv", aspect="auto");plt.colorbar();plt.show()
    # jsq_plasma
    plt.matshow(np.log10(jsq_m*boost)[::,:125].T[::-1, ::], cmap="hsv", aspect="auto");plt.colorbar();plt.show()


def plot_psi_interpolated(mhd_code, extent=(0, 5, 0, 0.25), txt_dir='/home/ilya/github/mhd_transfer/Release',
                          across_fraction=0.5):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib import colors

    def plot_single(value, extent, label, savename, aroundzero=False, vmin=None, vmax=None, show=False):
        fig, axes = plt.subplots(1, 1)

        if aroundzero:
            if vmin is None and vmax is None:
                vmin = np.min(value)
                vmax = np.max(value)
            print("Vmax = ", vmax)
            print("Vmin = ", vmin)
            assert vmin < 0
            assert vmax > 0
            norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
            im = axes.matshow(value[::,:].T[::-1, ::], cmap="seismic", vmin=vmin, vmax=vmax, norm=norm, aspect="auto", extent=extent)
        else:
            im = axes.matshow(value[::,:].T[::-1, ::], cmap="gist_rainbow", aspect="auto", vmin=vmin, vmax=vmax, extent=extent)

        axes.set_xlabel("z, pc")
        axes.set_ylabel("r, pc")
        divider = make_axes_locatable(axes)
        cax = divider.append_axes("right", size="5%", pad=0.00)
        cb = fig.colorbar(im, cax=cax)
        cb.set_label(label)
        fig.savefig(savename, bbox_inches="tight", dpi=300)
        if show:
            plt.show()

    Psi = np.loadtxt(os.path.join(txt_dir, "{}_Psi_psi_interpolated.txt".format(mhd_code)))
    shape = Psi.shape
    across = int(across_fraction*shape[1])
    Psi_m = np.ma.array(Psi, mask = Psi == 1)
    plot_single(Psi_m[:, :across], extent=extent, label=r"$\Psi$", savename="{}_psi_psi_interpolation.png".format(mhd_code))

    sigma = np.loadtxt(os.path.join(txt_dir, "{}_sigma_psi_interpolated.txt".format(mhd_code)))
    sigma_m = np.ma.array(sigma, mask = Psi == 1)
    plot_single(sigma_m[:, :across], extent=extent, label=r"$\sigma$", savename="{}_sigma_psi_interpolation.png".format(mhd_code))

    angle = np.loadtxt(os.path.join(txt_dir, "{}_angle_interpolated.txt".format(mhd_code)))
    angle_m = np.ma.array(angle, mask = Psi == 1)
    plot_single(np.rad2deg(angle_m[:, :across]), extent=extent, label=r"$\theta_{\rm pol}$, degree", vmax=8, savename="{}_poloidal_angle_interpolated.png".format(mhd_code))

    # theta = np.loadtxt(os.path.join(txt_dir, "{}_theta_psi_interpolated.txt".format(mhd_code)))
    # theta_m = np.ma.array(theta, mask = Psi == 1)
    # plot_single(np.rad2deg(theta_m), extent=extent, label=r"$\theta_{\rm pol}$, degree", vmax=15, savename="{}_poloidal_angle_interpolation.png".format(mhd_code))

    jsq_phi = np.loadtxt(os.path.join(txt_dir, "{}_jsq_phi_psi_interpolated.txt".format(mhd_code)))
    jsq_phi_m = np.ma.array(jsq_phi, mask=jsq_phi == 0)
    plot_single(np.log10(jsq_phi_m[:, :across]), extent=extent, label=r"$\lg{j_{\phi}}^2$, rel.units", savename="{}_jsq_phi_psi_interpolation.png".format(mhd_code))

    jsq = np.loadtxt(os.path.join(txt_dir, "{}_jsq_psi_interpolated.txt".format(mhd_code)))
    jsq_m = np.ma.array(jsq, mask=jsq == 0)
    plot_single(np.log10(jsq_m[:, :across]), extent=extent, label=r"$\lg{j_{\rm plasma}}^2$, rel.units", savename="{}_jsq_psi_interpolation.png".format(mhd_code))

    bphi = np.loadtxt(os.path.join(txt_dir, "{}_beta_phi_psi_interpolated.txt".format(mhd_code)))
    bphi_m = np.ma.array(bphi, mask=Psi == 1)
    # FIXME: in Lena's data beta_phi is nonnegative...
    plot_single(bphi_m[:, :across], extent=extent, label=r"$\beta_{\phi}$, c", savename="{}_beta_phi_psi_interpolation.png".format(mhd_code), aroundzero=False)
    # plot_single(bphi_m, extent=extent, label=r"$\beta_{\phi}$, c", savename="beta_phi_psi_interpolation.png",
    #             aroundzero=True, vmin=-0.1, vmax=0.1)

    Gamma = np.loadtxt(os.path.join(txt_dir, "{}_Gamma_psi_interpolated.txt".format(mhd_code)))
    Gamma_m = np.ma.array(Gamma, mask = Gamma == 1)
    plot_single(Gamma_m[:, :across], extent=extent, label=r"$\Gamma$", savename="{}_Gamma_psi_interpolation.png".format(mhd_code))

    N = np.loadtxt(os.path.join(txt_dir, "{}_N_psi_interpolated.txt".format(mhd_code)))
    N_m = np.ma.array(N, mask = N == 0)
    plot_single(np.log10(N_m[:, :across]), extent=extent, label=r"$\lg{N_{\rm plasma}}$, cm$^{-3}$", savename="{}_N_plasma_psi_interpolation.png".format(mhd_code))

    # FIXME: Lena's B_phi is negative
    B_phi = np.loadtxt(os.path.join(txt_dir, "{}_B_phi_psi_interpolated.txt".format(mhd_code)))
    B_phi_m = np.ma.array(B_phi, mask = Psi == 1)
    plot_single(np.log10(abs(B_phi_m[:, :across])), extent=extent, label=r"$\lg{B_{\phi}}$, G", savename="{}_B_phi_psi_interpolation.png".format(mhd_code))

    B_p = np.loadtxt(os.path.join(txt_dir, "{}_B_p_psi_interpolated.txt".format(mhd_code)))
    B_p_m = np.ma.array(B_p, mask = B_p == 0)
    plot_single(np.log10(B_p_m[:, :across]), extent=extent, label=r"$\lg{B_{\rm p}}$, G", savename="{}_B_p_psi_interpolation.png".format(mhd_code))

    B_tot = np.hypot(B_p_m, B_phi_m/Gamma_m)
    plot_single(np.log10(B_tot[:, :across]), extent=extent, label=r"$\lg{B_{\rm tot,plasma}}$, G", savename="{}_B_tot_plasma_psi_interpolation.png".format(mhd_code))


    # plt.matshow(np.log((B_tot**2)[::,:125].T[::-1, ::]), cmap="hsv", aspect="auto");plt.colorbar();plt.show()

    beta = np.sqrt(Gamma_m**2-1)/Gamma_m
    D = 1/(Gamma_m*(1-beta*np.cos(np.deg2rad(17.) + angle_m)))
    boost = D**2.5
    # plt.matshow(boost[::,:125].T[::-1, ::], cmap="hsv", aspect="auto");plt.colorbar();plt.show()
    plot_single(boost[:, :across], extent=extent, label=r"Boosting factor", savename="{}_boost_psi_interpolation.png".format(mhd_code))

    # n_plasma
    # plt.matshow(np.log10(N_m*boost/Gamma_m)[::,:125].T[::-1, ::], cmap="hsv", aspect="auto");plt.colorbar();plt.show()
    # B_tot
    # plt.matshow(np.log10(B_tot*boost)[::,:125].T[::-1, ::], cmap="hsv", aspect="auto");plt.colorbar();plt.show()
    # jsq_plasma
    # plt.matshow(np.log10(jsq_m*boost)[::,:125].T[::-1, ::], cmap="hsv", aspect="auto");plt.colorbar();plt.show()


def plot_images(mhd_code, rt_code, txt_files_dir="/home/ilya/github/mhd_transfer/Release", save_dir=None,
                n_along=600, n_across=200, lg_pixel_size_mas_min=np.log10(0.01), lg_pixel_size_mas_max=np.log10(0.1),
                freq_ghz_high=15.4, freq_ghz_low=None, cmap="magma", plot_beta_app=False, show=True):
    if save_dir is None:
        save_dir = os.getcwd()

    if freq_ghz_low:
        i_image_low = np.loadtxt("{}/{}_jet_image_{}_{}.txt".format(txt_files_dir, mhd_code, "i", freq_ghz_low))
    i_image_high = np.loadtxt("{}/{}_jet_image_{}_{}.txt".format(txt_files_dir, mhd_code, "i", freq_ghz_high))

    if freq_ghz_low:
        alpha_image = np.log(i_image_low/i_image_high)/np.log(freq_ghz_low/freq_ghz_high)

    jm = JetImage(z=0.00436, n_along=n_along, n_across=n_across,
                  lg_pixel_size_mas_min=lg_pixel_size_mas_min,
                  lg_pixel_size_mas_max=lg_pixel_size_mas_max, jet_side=True)

    jm.load_image_stokes("I", "{}/{}_jet_image_i_{}.txt".format(txt_files_dir, mhd_code, freq_ghz_high))
    if freq_ghz_low:
        jm.load_image_alpha(alpha_image)
    fig = jm.plot(log=True, Nan2zero=False, zoom_fr=1.0, axis_units="mas", figsize=(20, 7.5), cmap=cmap)
    fig.savefig(os.path.join(save_dir, "I_freq_{}_GHz_mhd_{}_rt_{}.png".format(freq_ghz_high, mhd_code, rt_code)),
                bbox_inches="tight", dpi=300)
    if show:
        plt.show()
    plt.close(fig)

    if plot_beta_app:
        jm.load_image_stokes("beta_app", "{}/{}_jet_image_beta_app_{}.txt".format(txt_files_dir, mhd_code, freq_ghz_high))
        fig = jm.plot(log=False, Nan2zero=False, zoom_fr=1.0, axis_units="mas", figsize=(20, 7.5), cmap="gist_ncar_r",
                      cblabel=r"$\beta_{\rm app}$, c", scale_by_pixsize=False, beta_app_min=0, beta_app_max=7)
        fig.savefig(os.path.join(save_dir, "beta_app_freq_{}_GHz_mhd_{}_rt_{}.png".format(freq_ghz_high, mhd_code, rt_code)),
                    bbox_inches="tight", dpi=300)
        if show:
            plt.show()
        plt.close(fig)

    if freq_ghz_low:
        fig = jm.plot_alpha(figsize=(20, 7.5), alpha_min=-1., alpha_max=0.0, count_levels_from_image_min=True)
        fig.savefig(os.path.join(save_dir, "alpha_mhd_{}_rt_{}.png".format(mhd_code, rt_code)),
                    bbox_inches="tight", dpi=300)
        if show:
            plt.show()
        plt.close(fig)


def plot_several_images(prefix, directory):
    import glob
    i_images_files = sorted(glob.glob(os.path.join(directory, prefix + "_*_jet_image_i_15.4.txt")))
    for i_image_file in i_images_files:
        i_image_file = os.path.split(i_image_file)[-1]
        print("Plotting ", i_image_file)
        psi = i_image_file.split("_")[2]
        dpsi = i_image_file.split("_")[4]
        subcode = "_psi_{}_dpsi_{}".format(psi, dpsi)
        code = prefix + subcode
        plot_images(code, "byhand", txt_files_dir=directory, save_dir=directory, n_along=600, n_across=100,
                    lg_pixel_size_mas_min=np.log10(0.01), lg_pixel_size_mas_max=np.log10(0.05), freq_ghz_high=15.4,
                    plot_beta_app=True, show=False)


def plot_raw(txt, label, extent, savename=None):
    toplot = np.loadtxt(txt)
    fig, axes = plt.subplots(1, 1)
    im = axes.matshow(toplot, extent=extent, aspect="equal")
    axes.set_xlabel(r"$z_{\rm obs}$, mas")
    axes.set_ylabel(r"$r$, mas")
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(axes)
    cax = divider.append_axes("right", size="5%", pad=0.00)
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(label)
    if savename is not None:
        fig.savefig(savename, bbox_inches="tight", dpi=300)
    plt.show()


if __name__ == "__main__":

    # z = 0.00436
    # freq_ghz = 15.4
    # number_of_pixels_along = 400
    # number_of_pixels_across = 100
    # pixel_size_mas_start = 0.05
    # pixel_size_mas_stop = 0.05
    # txt_dir = "/home/ilya/github/mhd_transfer/cmake-build-debug"
    # save_dir = "/home/ilya/data/Lena/RSF2022"
    #
    # model_stokes_dict = {"Q": os.path.join(txt_dir, "m1s200g2b881.675162r0.000312_n_none_jet_image_q_15.4.txt"),
    #                      "U": os.path.join(txt_dir, "m1s200g2b881.675162r0.000312_n_none_jet_image_u_15.4.txt")}
    #
    # jm = JetImage(z,
    #               number_of_pixels_along, number_of_pixels_across,
    #               np.log10(pixel_size_mas_start), np.log10(pixel_size_mas_stop))
    # jm.load_image_stokes("I", os.path.join(txt_dir, "m1s200g2b881.675162r0.000312_n_none_jet_image_i_15.4.txt"))
    # jm.plot(save_picture_file_name="test.png", save_dir=save_dir, figsize=(10, 5), plot_polarization=True,
    #         model_txt_files_dict=model_stokes_dict)


    # This applies to my local work in CLion. Change it to ``Release`` (or whatever) if necessary.
    jetpol_run_directory = "/home/ilya/github/mhd_transfer/cmake-build-debug"
    n_along = 400
    n_across = 100
    lg_pixel_size_mas_min = np.log10(0.05)
    # lg_pixel_size_mas_min = np.log10(0.003)
    # lg_pixel_size_mas_max = np.log10(0.1)
    lg_pixel_size_mas_max = np.log10(0.05)
    freq_ghz_high = 15.4
    freq_ghz_low = None
    # mhd_code = "best"
    # mhd_code = "m2s10g2b44.614955r0.000595"
    # mhd_code = "m1s10g2b123.971372r0.000369"
    mhd_code = "m1s200g2b881.675162r0.000312_n_none"
    # mhd_code = "m1s10g2b123.971372r0.000369_psi_0.750000_dpsi_0.015000"
    rt_code = "n"

    plot_images(mhd_code=mhd_code, rt_code=rt_code, freq_ghz_high=freq_ghz_high, freq_ghz_low=freq_ghz_low,
                n_along=n_along, n_across=n_across, save_dir="/home/ilya/data/Lena/RSF2022",
                lg_pixel_size_mas_min=lg_pixel_size_mas_min, lg_pixel_size_mas_max=lg_pixel_size_mas_max,
                cmap="jet", plot_beta_app=False, txt_files_dir=jetpol_run_directory)

    import sys; sys.exit(0)


    # # Test case - just plotting picture
    # # stokes = ("I", "Q", "U", "V")
    # stokes = ("I",)
    # # FIXME: Substitute with values used in radiative transfer
    # cjms = [JetImage(z=0.00436, n_along=200, n_across=50, lg_pixel_size_mas_min=np.log10(0.05),
    #                  lg_pixel_size_mas_max=np.log10(0.05), jet_side=False) for _ in stokes]
    # for i, stk in enumerate(stokes):
    #     cjms[i].load_image_stokes(stk, "{}/cjet_image_{}.txt".format(jetpol_run_directory, stk.lower()))
    #
    # jms = [JetImage(z=0.00436, n_along=200, n_across=50, lg_pixel_size_mas_min=np.log10(0.05),
    #                 lg_pixel_size_mas_max=np.log10(0.05), jet_side=True) for _ in stokes]
    # for i, stk in enumerate(stokes):
    #     jms[i].load_image_stokes(stk, "{}/jet_image_{}.txt".format(jetpol_run_directory, stk.lower()))
    #
    # j = TwinJetImage(jms[0], cjms[0])
    # j.plot_contours(zoom_fr=0.2, nlevels=20, aspect="auto")
    # plt.show()

    # # Test case - convolving model image with beam (for comparing with CLEAN maps)
    # stokes = ("I", "Q", "U")
    # jms = [JetImage(z=0.05, n_along=600, n_across=200, lg_pixel_size_mas_min=-3, lg_pixel_size_mas_max=-1, jet_side=True) for _ in stokes]
    # for i, stk in enumerate(stokes):
    #     jms[i].load_image_stokes(stk, "../{}/jet_image_{}.txt".format(jetpol_run_directory, stk.lower()))
    #     jms[i].save_image_to_difmap_format("dfm_{}.mod".format(stk))
    #     convert_difmap_model_file_to_CCFITS("dfm_{}.mod".format(stk), stk, (512, 0.1), (0.6, 0.6, 0), "1458+718.u.2006_09_06.uvf",
    #                                         "convolved_{}.fits".format(stk))

