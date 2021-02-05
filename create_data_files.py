import os
import glob
import numpy as np
import matplotlib.pyplot as plt

files_dir = "/home/ilya/github/mhd_solver/Release"

r_pc_files = glob.glob(os.path.join(files_dir, "*_r_pc.txt"))
n_plasma_files = glob.glob(os.path.join(files_dir, "*_n_plasma.txt"))
Gamma_files = glob.glob(os.path.join(files_dir, "*_Gamma.txt"))
Bphi_files = glob.glob(os.path.join(files_dir, "*_B_phi.txt"))
Bp_files = glob.glob(os.path.join(files_dir, "*_B_p.txt"))
beta_phi_files = glob.glob(os.path.join(files_dir, "*_beta_phi.txt"))
I_files = glob.glob(os.path.join(files_dir, "*_I.txt"))

assert len(I_files) == len(beta_phi_files) == len(Bp_files) == len(Bphi_files) == len(Gamma_files) == len(n_plasma_files) == len(r_pc_files)
n_files = len(r_pc_files)

# Find distance from BH where profiles are calculated
zs_pc = sorted([float(os.path.split(path)[-1].split("_")[0]) for path in r_pc_files])

result_z_pc = list()
result_r_pc = list()
result_n_plasma = list()
result_Gamma = list()
result_Bphi = list()
result_Bp = list()
result_beta_phi = list()
result_jsq_plasma = list()

fig, axes = plt.subplots(1, 1)
for i, z_pc in enumerate(zs_pc):
    r_pc = np.loadtxt(os.path.join(files_dir, "{:.6f}_r_pc.txt".format(z_pc)))
    n_plasma = np.loadtxt(os.path.join(files_dir, "{:.6f}_n_plasma.txt".format(z_pc)))
    Gamma = np.loadtxt(os.path.join(files_dir, "{:.6f}_Gamma.txt".format(z_pc)))
    Bphi = np.loadtxt(os.path.join(files_dir, "{:.6f}_B_phi.txt".format(z_pc)))
    Bp = np.loadtxt(os.path.join(files_dir, "{:.6f}_B_p.txt".format(z_pc)))
    beta_phi = np.loadtxt(os.path.join(files_dir, "{:.6f}_beta_phi.txt".format(z_pc)))
    I = np.loadtxt(os.path.join(files_dir, "{:.6f}_I.txt".format(z_pc)))
    jsq_plasma = (np.gradient(I, r_pc, edge_order=2)/r_pc)**2 * Gamma**2
    axes.plot(r_pc, jsq_plasma, '.', label="z = {:.3f} pc".format(z_pc))

    result_z_pc.extend(len(r_pc)*[z_pc])
    result_r_pc.extend(r_pc)
    result_n_plasma.extend(n_plasma)
    result_Gamma.extend(Gamma)
    result_Bphi.extend(Bphi)
    result_Bp.extend(Bp)
    result_beta_phi.extend(beta_phi)
    result_jsq_plasma.extend(jsq_plasma)

axes.set_xlabel("R, pc")
# plt.yscale("log")
# plt.xscale("log")
plt.legend()
plt.show()

result_jsq_plasma = np.array(result_jsq_plasma)
result_jsq_plasma /= np.max(result_jsq_plasma)


np.savetxt(os.path.join(files_dir, "Gamma_field.txt"), np.dstack((result_z_pc, result_r_pc, result_Gamma))[0])
np.savetxt(os.path.join(files_dir, "n_plasma_field.txt"), np.dstack((result_z_pc, result_r_pc, result_n_plasma))[0])
np.savetxt(os.path.join(files_dir, "B_p_field.txt"), np.dstack((result_z_pc, result_r_pc, result_Bp))[0])
np.savetxt(os.path.join(files_dir, "B_phi_field.txt"), np.dstack((result_z_pc, result_r_pc, result_Bphi))[0])
np.savetxt(os.path.join(files_dir, "beta_phi_field.txt"), np.dstack((result_z_pc, result_r_pc, result_beta_phi))[0])
np.savetxt(os.path.join(files_dir, "jsq_plasma_field.txt"), np.dstack((result_z_pc, result_r_pc, result_jsq_plasma))[0])