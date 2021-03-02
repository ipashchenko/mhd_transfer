import sys
import os
import numpy as np
from create_data_files_psi import create_files_for_transfer
from jet_image import plot_images
import dlib


# Speed of light [cm / s]
c = 29979245800.0
# Gravitational constant [cm3 / (g s2)]
G = 6.67430e-08
# Mass of the Sun [g]
M_sun = 1.98840987e+33

M_BH = 6.5E+09
Psi_tot = 1E+34
r_g = G*M_BH/c**2


solver_dir = "/home/ilya/github/mhd_solver/Release"
transfer_dir = "/home/ilya/github/mhd_transfer/Release"
images_save_dir = None
executable_sf = None
executable_mf = None


def find_scale(mhd_code, executable_sf, particles_heating_model, required_flux_Jy=1.0, gamma_min=100.0, executable_dir=None, changing_s=False):
    cwd = os.getcwd()
    if executable_dir is None:
        executable_dir = cwd
        os.chdir(executable_dir)

    if changing_s:
        changing_s = "true"
    else:
        changing_s = "false"

    def flux_diff(scale):
        os.system("./{} {} {} {} {} {}".format(executable_sf, mhd_code, scale, gamma_min, changing_s, particles_heating_model))
        image = np.loadtxt(mhd_code + "_jet_image_i_15.4.txt")
        abs_log_diff = abs(np.log(required_flux_Jy) - np.log(np.sum(image)))
        print("==========================================================")
        print("For scale = {:.3f}, lg(Flux) = {:.3f}, diff log = {:.3f}".format(scale, np.log10(np.sum(image)), abs_log_diff))
        print("==========================================================")
        if np.isnan(abs_log_diff):
            abs_log_diff = 1e3
        return abs_log_diff

    lower_bounds = [1e-03]
    upper_bounds = [1e+10]
    n_fun_eval = 10
    result = dlib.find_min_global(flux_diff, lower_bounds, upper_bounds, n_fun_eval)
    os.chdir(cwd)

    return result


for sigma_M in np.logspace(np.log10(5), np.log10(50), 5, base=10):
    print("sigma_M = {:.1f}".format(sigma_M))
    for B_L in np.logspace(np.log10(1), np.log10(100), 5, base=10):
        print("B_L = {:.1f}".format(B_L))
        R_L = np.sqrt(Psi_tot/(np.pi*B_L))
        r_rato = r_g/R_L
        a = 8*r_rato/(1 + 16*r_rato**2)
        for gamma_in in (2, 3.4, 5):
            print("gamma_in = {:.1f}".format(gamma_in))
            mhd_code = "m1s{}g{}b{}".format(int(sigma_M), int(gamma_in), int(B_L))
            os.chdir(solver_dir)
            # Solve MHD equations
            os.system("./mhd_solver {} {} {} {} {}".format(mhd_code, sigma_M, B_L, a, gamma_in))

            # Re-format profiles files to suitable for RT
            create_files_for_transfer(mhd_code, solver_dir, transfer_dir)

            # Radiation Transport
            os.chdir(transfer_dir)
            for rt in ("n", "bsq", "jsq"):
                print("Particles heating model = {}".format(rt))
                for gamma_min in (1, 100, 200):
                    print("gamma_min = {:.1f}".format(gamma_min))
                    rt_code = "{}_{}"
                    # Find scale factor of NT particles
                    res = find_scale(mhd_code, executable_sf, particles_heating_model=rt, required_flux_Jy=1.0,
                                     gamma_min=gamma_min, executable_dir=transfer_dir, changing_s=False)
                    scale = res[0][0]
                    os.system("./{} {} {} {} {} {}".format(executable_mf, mhd_code, scale, gamma_min, "false", rt))
                    plot_images(mhd_code=mhd_code, rt_code=rt_code, txt_files_dir=transfer_dir, save_dir=images_save_dir)

                    sys.exit(0)