from __future__ import print_function
import matplotlib
matplotlib.use('Agg')

import os
import sys

import argparse

import json
import numpy as np

from ase.neb import fit0
from ase.io import write
from ase.optimize import FIRE
from ase.optimize.precon import PreconFIRE, Exp, PreconLBFGS


from distutils import spawn

from imeall import app
from ForceMixerCarver import ForceMixingCarvingCalculator
from matscipy.socketcalc import VaspClient, SocketCalculator

from quippy import AtomsReader, AtomsWriter, Atoms
from quippy import Potential

import matplotlib.pyplot as plt

from imeall.calc_elast_dipole import find_h_atom

#from ase.neb import NEB
from nebForceIntegrator import NEB

class NEBAnalysis(object):
    def __init__(self):
        pass

    def save_barriers(self, images, neb, prefix=""):
        """ Saves the NEB results and plots energy barrier.
        Arguments:
          images ::
          prefix ::
          neb :: :class:ase:`NEB`
        """
        if not os.path.exists('NEB_images'):
            os.mkdir('NEB_images')

        for i, image in enumerate(images):
            write("NEB_images/" + prefix + "ini_NEB_image%03i.xyz" % i, image)

        plot_dir = "NEB_plots"
        if not os.path.exists(plot_dir):
            os.mkdir(plot_dir)

        fig = plt.figure()
        ax = fig.add_subplot(111)

        R = [atoms.positions for atoms in images]
        E = neb.get_potential_energies()
        F = [atoms.get_forces() for atoms in images]
        A = images[0].cell
        pbc = images[0].pbc
        s, E, Sfit, Efit, lines = fit0(E, F, R, A, pbc)

        s = np.array(s)
        norm = s.max()
        s /= norm
        Sfit /= norm

        for x, y in lines:
            x /= norm
            ax.plot(x, y, '-C0')

        ax.plot(Sfit, Efit, 'C0-')
        Ef = max(Efit) - E[0]
        ax.plot(s, E, "oC0", label="%s: $E_f = %.4f$ (meV)"%("QM/MM EAM result", Ef*1000))
        ax.legend(loc="best")
        ax.grid(True, linestyle="dashed")

        np.savetxt(plot_dir + '/' + prefix + "_s_E.txt", np.array([s, E]))
        np.savetxt(plot_dir + '/' + prefix + "_fit.txt", np.array([Sfit, Efit]))
        np.savetxt(plot_dir + '/' + prefix + "_lines.txt", np.vstack(lines))

        fig.savefig(plot_dir + '/' + prefix + "_barrier.eps")
        return None


class NEBPaths(object):
    def __init__(self):
        pass

    def build_h_nebpath(self, struct_path="fe_bcc.xyz", neb_path=np.array([0.25, 0.0, -0.25]), alat=2.8297, knots=5, sup_cell=5, fmax=1.e-4):
        """
        Takes a vector neb_path, and generates n intermediate images along the minimum energy path.
        the struct path should point to the relaxed structure.
        """
        POT_DIR = os.path.join(app.root_path, 'potentials')
        eam_pot = os.path.join(POT_DIR, 'PotBH.xml')
        r_scale = 1.00894848312
        mm_pot = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)

        ats_ini = AtomsReader(struct_path)[-1]
        #Pick the top tetrahedral H position in the cell [[1,0,0],[0,1,0],[0,0,1]]
        tetra_pos = alat*np.array([0.5, 0.0, 0.75])
        mid_point = 0.5*(np.diag(ats_ini.get_cell()))
        mid_point = [((sup_cell-1)/2.)*alat for sp in range(3)]
        h_pos = tetra_pos + mid_point

        disloc_ini = ats_ini.copy()
        h_tmp = h_pos
        disloc_ini.add_atoms(np.array(h_tmp), 1)

        disloc_fin = ats_ini.copy()
        h_tmp = h_pos + neb_path*alat
        disloc_fin.add_atoms(h_tmp,1)

#Relax images at the start and end of trajectory.
        disloc_ini.set_calculator(mm_pot)
        opt = FIRE(disloc_ini)
        opt.run(fmax=fmax)

        disloc_fin.set_calculator(mm_pot)
        opt = FIRE(disloc_fin)
        opt.run(fmax=fmax)

        return disloc_ini, disloc_fin

if __name__=="__main__":
    mpirun = spawn.find_executable('mpirun')
    vasp = '/home/mmm0007/vasp/vasp.5.4.1/bin/vasp_std'

    parser = argparse.ArgumentParser()
    parser.add_argument("--buff", "-b", type=float, default=8.0)
    parser.add_argument("--qm_radius", "-q", type=float, default=2.0)
    parser.add_argument("--use_socket", "-u", action="store_false")
    parser.add_argument("--neb_path", "-n", nargs="+", type=float, default=[0.25, 0, -0.25])
    parser.add_argument("--fmax", "-f", type=float, help="maximum force for relaxation.", default=0.08)
    parser.add_argument("--sup_cell", "-s", type=int, help="size of fe matrix super cell")
    parser.add_argument("--knots", "-kn", type=int, help="number of images", default=15)
    parser.add_argument("--k", "-k", type=float, default=5.0, help="spring constant for NEB (default eV/A^2).")
    parser.add_argument("--input_file", "-i", default="fe_bcc_h.xyz", help="input file.")
    parser.add_argument("--auto_gen", "-a", action="store_true")
    parser.add_argument("--use_gap", "-g", action="store_true")
    args = parser.parse_args()

    #only one qmpot specifed at command line
    use_eampot = not(args.use_gap or args.use_socket) 
    assert sum(args.use_gap + args.use_socket + use_mmpot) == 1 

    buff = args.buff
    qm_radius = args.qm_radius

    POT_DIR = os.path.join(app.root_path, 'potentials')
    eam_pot = os.path.join(POT_DIR, 'PotBH.xml')
    r_scale = 1.00894848312
    mm_pot = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)

    if args.auto_gen:
        gb_cell = AtomsReader(args.input_file)[-1]
    else:
        gb_cell = AtomsReader("disloc_ini_traj.xyz")[-1]

    defect = find_h_atom(gb_cell)
    h_pos = defect.position
    x, y, z = gb_cell.positions.T
    radius1 = np.sqrt((x - h_pos[0])**2 + (y-h_pos[1])**2 + (z-h_pos[2])**2)

    qm_region_mask = (radius1 < qm_radius)
    qm_buffer_mask = (radius1 < qm_radius + buff)

    print ("\nNumber of atoms in qm region of %.1f" % qm_radius +
                                    "A : %i" % np.count_nonzero(qm_region_mask))

    print ("together with the buffer of %.1f" % (qm_radius + buff ) +
                                    "A %i" % np.count_nonzero(qm_buffer_mask))

    magmoms=[2.6, len(gb_cell)]
    vasp_args = dict(xc='PBE', amix=0.01, amin=0.001, bmix=0.001, amix_mag=0.01, bmix_mag=0.001,
                     kpts=[1, 1, 1], kpar=1, lreal='auto', nelmdl=-15, ispin=2, prec='Accurate', ediff=1.e-4,
                     nelm=100, algo='VeryFast', lplane=False, lwave=False, lcharg=False, istart=0, encut=400,
                     magmom=magmoms, maxmix=30, voskown=0, ismear=1, sigma=0.1, isym=0) # possibly try iwavpr=12, should be faster if it works

    # parallel config.
    procs = 48
    kpts = [1, 1, 1]
    # need to have procs % n_par == 0
    n_par = 1
    if procs <= 8:
        n_par = procs
    else:
        for _npar in range(2, int(np.sqrt(1.*procs))):
            if procs % int(_npar) == 0:
                n_par = procs // int(_npar)

    if args.use_socket:
        vasp_client = VaspClient(client_id=0, npj=procs, ppn=1,
                                 exe=vasp, mpirun=mpirun, parmode='mpi',
                                 ibrion=13, nsw=1000000,
                                 npar=n_par, **vasp_args)
        qm_pot = SocketCalculator(vasp_client)
    elif args.use_gap:
        gap_pot = os.path.join(POT_DIR,'gp33b.xml')
        sparse_file = 'gp33b.xml.sparseX.GAP_2016_10_3_60_19_29_10_8911'
        gap_pot_sparse = os.path.join(POT_DIR, sparse_file)
        shutil.copy(gap_pot, './')
        shutil.copy(gap_pot_sparse, './')
        qm_pot = Potential('IP GAP', param_filename=gap_pot)
    else:
    #for entirely mm potential
        qm_pot = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(1.00894848312), param_filename=eam_pot)

    nebpath = NEBPaths()
    nebanalysis = NEBAnalysis()

    if args.auto_gen:
        disloc_ini, disloc_fin = nebpath.build_h_nebpath(neb_path=np.array(args.neb_path), fmax = args.fmax, sup_cell=args.sup_cell)
    else:
        disloc_ini = Atoms('disloc_ini_traj.xyz')
        disloc_fin = Atoms('disloc_fin_traj.xyz')

    n_knots = args.knots
    images = [disloc_ini] + \
             [disloc_ini.copy() for i in range(n_knots)] + \
             [disloc_fin]

    qmmm_pot = ForceMixingCarvingCalculator(disloc_ini, qm_region_mask,
                                            mm_pot, qm_pot,
                                            buffer_width=buff,
                                            pbc_type=[False, False, False])

    for image in images:
        image.set_calculator(qmmm_pot)

    neb = NEB(images, k=args.k, force_only=True)
    neb.interpolate(mic=False)
    #nebanalysis.save_barriers(images, neb, prefix="ini")
    #opt = FIRE(neb)
    opt = PreconLBFGS(neb)
    opt.run(fmax=args.fmax)

    nebanalysis.save_barriers(images, neb, prefix="fin")

    if args.use_socket:
        sock_calc.shutdown()
