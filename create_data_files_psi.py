# This creates data files in format (z[pc], Psi, value) and single (z[pc], r[pc], Psi)
import os
import sys
import glob
import numpy as np
import matplotlib
label_size = 20
matplotlib.rcParams['xtick.labelsize'] = label_size
matplotlib.rcParams['ytick.labelsize'] = label_size
matplotlib.rcParams['axes.titlesize'] = label_size
matplotlib.rcParams['axes.labelsize'] = label_size
matplotlib.rcParams['font.size'] = label_size
matplotlib.rcParams['legend.fontsize'] = label_size
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def create_files_for_transfer(run_name, files_dir, save_dir):

    r_pc_files = glob.glob(os.path.join(files_dir, "{}_*_r_pc.txt".format(run_name)))
    n_plasma_files = glob.glob(os.path.join(files_dir, "{}_*_n_plasma.txt".format(run_name)))
    Gamma_files = glob.glob(os.path.join(files_dir, "{}_*_Gamma.txt".format(run_name)))
    Bphi_files = glob.glob(os.path.join(files_dir, "{}_*_B_phi.txt".format(run_name)))
    Bp_files = glob.glob(os.path.join(files_dir, "{}_*_B_p.txt".format(run_name)))
    Bsq_plasma_files = glob.glob(os.path.join(files_dir, "{}_*_Bsq_plasma.txt".format(run_name)))
    beta_phi_files = glob.glob(os.path.join(files_dir, "{}_*_beta_phi.txt".format(run_name)))
    I_files = glob.glob(os.path.join(files_dir, "{}_*_I.txt".format(run_name)))
    Psi_files = glob.glob(os.path.join(files_dir, "{}_*_Psi.txt".format(run_name)))

    assert len(I_files) == len(Psi_files) == len(beta_phi_files) == len(Bp_files) == len(Bphi_files) == len(Bsq_plasma_files) == len(Gamma_files) == len(n_plasma_files) == len(r_pc_files)
    n_files = len(r_pc_files)

    # Find distance from BH where profiles are calculated
    zs_pc = sorted([float(os.path.split(path)[-1].split("_")[1]) for path in r_pc_files])

    result_z_pc = list()
    result_r_pc = list()
    result_n_plasma = list()
    result_Gamma = list()
    result_Bphi = list()
    result_Bp = list()
    result_Bsq_plasma = list()
    result_Psi = list()
    result_beta_phi = list()
    result_jsq_plasma = list()

    for i, z_pc in enumerate(zs_pc):
        r_pc = np.loadtxt(os.path.join(files_dir, "{}_{:.6f}_r_pc.txt".format(run_name, z_pc)))
        n_plasma = np.loadtxt(os.path.join(files_dir, "{}_{:.6f}_n_plasma.txt".format(run_name, z_pc)))
        Gamma = np.loadtxt(os.path.join(files_dir, "{}_{:.6f}_Gamma.txt".format(run_name, z_pc)))
        Bphi = np.loadtxt(os.path.join(files_dir, "{}_{:.6f}_B_phi.txt".format(run_name, z_pc)))
        Bp = np.loadtxt(os.path.join(files_dir, "{}_{:.6f}_B_p.txt".format(run_name, z_pc)))
        Bsq_plasma = np.loadtxt(os.path.join(files_dir, "{}_{:.6f}_Bsq_plasma.txt".format(run_name, z_pc)))
        beta_phi = np.loadtxt(os.path.join(files_dir, "{}_{:.6f}_beta_phi.txt".format(run_name, z_pc)))
        I = np.loadtxt(os.path.join(files_dir, "{}_{:.6f}_I.txt".format(run_name, z_pc)))
        Psi = np.loadtxt(os.path.join(files_dir, "{}_{:.6f}_Psi.txt".format(run_name, z_pc)))
        jsq_plasma = (np.gradient(I, r_pc, edge_order=2)/r_pc)**2 * Gamma**2

        result_z_pc.extend(len(r_pc)*[z_pc])
        result_r_pc.extend(r_pc)
        result_n_plasma.extend(n_plasma)
        result_Gamma.extend(Gamma)
        result_Bphi.extend(Bphi)
        result_Bp.extend(Bp)
        result_Bsq_plasma.extend(Bsq_plasma)
        result_beta_phi.extend(beta_phi)
        result_jsq_plasma.extend(jsq_plasma)
        result_Psi.extend(Psi)

    result_jsq_plasma = np.array(result_jsq_plasma)
    result_jsq_plasma /= np.max(result_jsq_plasma)
    result_z_pc = np.array(result_z_pc)
    result_r_pc = np.array(result_r_pc)

    np.savetxt(os.path.join(save_dir, "{}_Gamma_field_psi.txt".format(run_name)), np.dstack((result_z_pc, result_Psi, result_Gamma))[0])
    np.savetxt(os.path.join(save_dir, "{}_n_plasma_field_psi.txt".format(run_name)), np.dstack((result_z_pc, result_Psi, result_n_plasma))[0])
    np.savetxt(os.path.join(save_dir, "{}_B_p_field_psi.txt".format(run_name)), np.dstack((result_z_pc, result_Psi, result_Bp))[0])
    np.savetxt(os.path.join(save_dir, "{}_Bsq_plasma_field_psi.txt".format(run_name)), np.dstack((result_z_pc, result_Psi, result_Bsq_plasma))[0])
    np.savetxt(os.path.join(save_dir, "{}_Psi_field.txt".format(run_name)), np.dstack((result_z_pc, result_r_pc, result_Psi))[0])
    np.savetxt(os.path.join(save_dir, "{}_Psi_field_psi.txt".format(run_name)), np.dstack((result_z_pc, result_Psi, result_Psi))[0])
    np.savetxt(os.path.join(save_dir, "{}_B_phi_field_psi.txt".format(run_name)), np.dstack((result_z_pc, result_Psi, result_Bphi))[0])
    np.savetxt(os.path.join(save_dir, "{}_beta_phi_field_psi.txt".format(run_name)), np.dstack((result_z_pc, result_Psi, result_beta_phi))[0])
    np.savetxt(os.path.join(save_dir, "{}_jsq_plasma_field_psi.txt".format(run_name)), np.dstack((result_z_pc, result_Psi, result_jsq_plasma))[0])


def create_files_for_transfer_Lena(run_name, files_dir, save_dir):

    r_pc_files = glob.glob(os.path.join(files_dir, "{}_*_r_pc.txt".format(run_name)))
    n_plasma_files = glob.glob(os.path.join(files_dir, "{}_*_n_plasma.txt".format(run_name)))
    Gamma_files = glob.glob(os.path.join(files_dir, "{}_*_Gamma.txt".format(run_name)))
    Bphi_files = glob.glob(os.path.join(files_dir, "{}_*_B_phi.txt".format(run_name)))
    Bp_files = glob.glob(os.path.join(files_dir, "{}_*_B_p.txt".format(run_name)))
    # Bsq_plasma_files = glob.glob(os.path.join(files_dir, "{}_*_Bsq_plasma.txt".format(run_name)))
    # beta_phi_files = glob.glob(os.path.join(files_dir, "{}_*_beta_phi.txt".format(run_name)))
    I_files = glob.glob(os.path.join(files_dir, "{}_*_jsq_plasma.txt".format(run_name)))
    Psi_files = glob.glob(os.path.join(files_dir, "{}_*_Psi.txt".format(run_name)))

    assert len(I_files) == len(Psi_files) == len(Bp_files) == len(Bphi_files) == len(Gamma_files) == len(n_plasma_files) == len(r_pc_files)
    n_files = len(r_pc_files)

    # Find distance from BH where profiles are calculated
    zs_pc = sorted([float(os.path.split(path)[-1].split("_")[1]) for path in r_pc_files])

    result_z_pc = list()
    result_r_pc = list()
    result_n_plasma = list()
    result_Gamma = list()
    result_Bphi = list()
    result_Bp = list()
    result_Bsq_plasma = list()
    result_Psi = list()
    result_beta_phi = list()
    result_jsq_plasma = list()

    for i, z_pc in enumerate(zs_pc):
        r_pc = np.loadtxt(os.path.join(files_dir, "{}_{:.6f}_r_pc.txt".format(run_name, z_pc)))
        n_plasma = np.loadtxt(os.path.join(files_dir, "{}_{:.6f}_n_plasma.txt".format(run_name, z_pc)))
        Gamma = np.loadtxt(os.path.join(files_dir, "{}_{:.6f}_Gamma.txt".format(run_name, z_pc)))
        Bphi = np.loadtxt(os.path.join(files_dir, "{}_{:.6f}_B_phi.txt".format(run_name, z_pc)))
        Bp = np.loadtxt(os.path.join(files_dir, "{}_{:.6f}_B_p.txt".format(run_name, z_pc)))
        # beta_phi = np.loadtxt(os.path.join(files_dir, "{}_{:.6f}_beta_phi.txt".format(run_name, z_pc)))
        jsq_plasma = np.loadtxt(os.path.join(files_dir, "{}_{:.6f}_jsq_plasma.txt".format(run_name, z_pc)))
        # Bsq_plasma = np.loadtxt(os.path.join(files_dir, "{}_{:.6f}_Bsq_plasma.txt".format(run_name, z_pc)))
        Psi = np.loadtxt(os.path.join(files_dir, "{}_{:.6f}_Psi.txt".format(run_name, z_pc)))

        Bsq_plasma = Bp**2 + (Bphi/Gamma)**2
        beta_phi = np.zeros(len(Gamma))

        result_z_pc.extend(len(r_pc)*[z_pc])
        result_r_pc.extend(r_pc)
        result_n_plasma.extend(n_plasma)
        result_Gamma.extend(Gamma)
        result_Bphi.extend(Bphi)
        result_Bp.extend(Bp)
        result_Bsq_plasma.extend(Bsq_plasma)
        result_beta_phi.extend(beta_phi)
        result_jsq_plasma.extend(jsq_plasma)
        result_Psi.extend(Psi)

    result_jsq_plasma = np.array(result_jsq_plasma)
    result_jsq_plasma /= np.max(result_jsq_plasma)
    result_z_pc = np.array(result_z_pc)
    result_r_pc = np.array(result_r_pc)

    np.savetxt(os.path.join(save_dir, "{}_Gamma_field_psi.txt".format(run_name)), np.dstack((result_z_pc, result_Psi, result_Gamma))[0])
    np.savetxt(os.path.join(save_dir, "{}_n_plasma_field_psi.txt".format(run_name)), np.dstack((result_z_pc, result_Psi, result_n_plasma))[0])
    np.savetxt(os.path.join(save_dir, "{}_B_p_field_psi.txt".format(run_name)), np.dstack((result_z_pc, result_Psi, result_Bp))[0])
    np.savetxt(os.path.join(save_dir, "{}_Bsq_plasma_field_psi.txt".format(run_name)), np.dstack((result_z_pc, result_Psi, result_Bsq_plasma))[0])
    np.savetxt(os.path.join(save_dir, "{}_Psi_field.txt".format(run_name)), np.dstack((result_z_pc, result_r_pc, result_Psi))[0])
    np.savetxt(os.path.join(save_dir, "{}_Psi_field_psi.txt".format(run_name)), np.dstack((result_z_pc, result_Psi, result_Psi))[0])
    np.savetxt(os.path.join(save_dir, "{}_B_phi_field_psi.txt".format(run_name)), np.dstack((result_z_pc, result_Psi, result_Bphi))[0])
    np.savetxt(os.path.join(save_dir, "{}_beta_phi_field_psi.txt".format(run_name)), np.dstack((result_z_pc, result_Psi, result_beta_phi))[0])
    np.savetxt(os.path.join(save_dir, "{}_jsq_plasma_field_psi.txt".format(run_name)), np.dstack((result_z_pc, result_Psi, result_jsq_plasma))[0])


if __name__ == "__main__":
    files_dir = "/home/ilya/github/mhd_solver/Release"
    save_dir = "/home/ilya/github/mhd_transfer/Release"
    if len(sys.argv) != 2:
        raise Exception("Provide MHD code as a single argument")
    run_name = sys.argv[1]
    create_files_for_transfer(run_name, files_dir, save_dir)
