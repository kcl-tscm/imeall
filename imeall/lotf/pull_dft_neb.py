import os
import sys
import glob
import numpy as np

from ase import io
from ase.neb import fit0
from ase.io.extxyz import write_xyz
from ase import units
from os_supp.pushd import pushd

import matplotlib.pyplot as plt


def calc_neb_barriers(prefix, images):
    R = [ats.positions for ats in images]
    E = [ats.get_potential_energy() for ats in images]
    F = [atoms.get_forces() for atoms in images]
    A = images[0].cell
    pbc = images[0].pbc
    s, E, Sfit, Efit, lines = fit0(E, F, R, A, pbc)

#Normalizes path length to 1.
    s = np.array(s)
    norm = s.max()
    s /= norm
    Sfit /= norm

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for x, y in lines:
        x /= norm
        ax.plot(x, y, '-C0')

    ax.plot(Sfit, Efit, 'C0-')
    Ef = max(Efit) - E[0]
    ax.plot(s, E, "oC0", label="%s: $E_f = %.4f$ (meV)"%("Interpolated Energy Barrier", Ef*1000))
    ax.legend(loc="best")
    ax.grid(True, linestyle="dashed")
    np.savetxt(prefix + "_s_E.txt", np.array([s, E]))
    np.savetxt(prefix + "_fit.txt", np.array([Sfit, Efit]))
    np.savetxt(prefix + "_lines.txt", np.vstack(lines))
    fig.savefig(prefix + "_barrier.eps")

if __name__=='__main__':
    prefix = sys.argv[1]
    job_dirs = glob.glob("{prefix}*".format(prefix=prefix))
    job_dirs = filter(os.path.isdir, job_dirs) 
    func = lambda x: int(x.split("_")[1])
    job_dirs.sort(key = func)
    
    images=[]
    print "Removing old neb path"
    if os.path.isfile("neb_path.xyz"):
        os.remove("neb_path.xyz")
    for job in job_dirs:
        print "pulling ", job
        with pushd(job) as ctx0:
            try:
                ats = io.read('OUTCAR', index='-1')
                ats.arrays['forces'] = ats.get_forces()
            except StopIteration:
                print "Could not pull OUTCAR for:", job
            else:
                images.append(ats)

    for image in images:
        write_xyz("neb_path.xyz",image, append=True)
    calc_neb_barriers(prefix,images)

