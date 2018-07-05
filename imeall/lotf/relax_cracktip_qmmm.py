from __future__ import print_function
import matplotlib
matplotlib.use('Agg')

import os
import sys
import shutil

import argparse

import json
import numpy as np

from ase.neb import fit0
from ase.io import write, Trajectory
from ase.optimize import FIRE
from ase.optimize.precon import PreconFIRE, Exp, PreconLBFGS
from ase.io.xyz import write_xyz


from distutils import spawn

from imeall import app
from ForceMixerCarver import ForceMixingCarvingCalculator
from matscipy.socketcalc import VaspClient, SocketCalculator

from quippy import AtomsReader, AtomsWriter, Atoms
from quippy import Potential, set_fortran_indexing

import matplotlib.pyplot as plt

from imeall.calc_elast_dipole import find_h_atom
from simulate_crack import fix_edges

#from ase.neb import NEB
from nebForceIntegrator import NEB

set_fortran_indexing(False)

if __name__=="__main__":
    mpirun = spawn.find_executable('mpirun')
    vasp = '/home/mmm0007/vasp/vasp.5.4.1/bin/vasp_std'

    parser = argparse.ArgumentParser()
    parser.add_argument("--buff", "-b", type=float, default=9.0)
    parser.add_argument("--qm_radius", "-q", type=float, default=3.0)
    parser.add_argument("--use_socket", "-u", action="store_true")
    parser.add_argument("--fmax", "-f", type=float, help="maximum force for relaxation.", default=0.1)
    parser.add_argument("--sup_cell", "-s", type=int, help="size of fe matrix super cell")
    parser.add_argument("--input_file", "-i", default="crack.xyz", help="input file.")
    parser.add_argument("--use_gap", "-g", action="store_true")
    parser.add_argument("--precon", "-p", action='store_true')
    args = parser.parse_args()

    #only one qmpot specifed at command line
    use_eampot = not (args.use_gap or args.use_socket)
    if use_eampot:
        print ('USING EAMPOT')

    buff = args.buff
    qm_radius = args.qm_radius

    POT_DIR = os.path.join(app.root_path, 'potentials')
    eam_pot = os.path.join(POT_DIR, 'PotBH.xml')
    r_scale = 1.00894848312
    mm_pot = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)

    crack_slab = AtomsReader('crack_sim.xyz')[-1]

    crack_pos = crack_slab.info['CrackPos']
    x, y, z = crack_slab.positions.T
    radius1 = np.sqrt((x - crack_pos[0])**2 + (y-crack_pos[1])**2 + (z-crack_pos[2])**2)

    qm_region_mask = (radius1 < qm_radius)
    qm_buffer_mask = (radius1 < qm_radius + buff)

    print ("\nNumber of atoms in qm region of %.1f" % qm_radius +
                                    "A : %i" % np.count_nonzero(qm_region_mask))

    print ("together with the buffer of %.1f" % (qm_radius + buff ) +
                                    "A %i" % np.count_nonzero(qm_buffer_mask))

    if args.use_socket:
        magmoms=[2.6, len(crack_slab)]
        vasp_args = dict(xc='PBE', amix=0.01, amin=0.001, bmix=0.001, amix_mag=0.01, bmix_mag=0.001,
                     kpts=[1, 1, 4], kpar=1, lreal='auto', nelmdl=-15, ispin=2, prec='Accurate', ediff=1.e-4,
                     nelm=100, algo='VeryFast', lplane=False, lwave=False, lcharg=False, istart=0, encut=400,
                     magmom=magmoms, maxmix=30, voskown=0, ismear=1, sigma=0.1, isym=0) # possibly try iwavpr=12, should be faster if it works

    # parallel config.
        procs = 96
    # need to have procs % n_par == 0
        n_par = 1
        if procs <= 8:
            n_par = procs
        else:
            for _npar in range(2, int(np.sqrt(1.*procs))):
                if procs % int(_npar) == 0:
                    n_par = procs // int(_npar)
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

    strain_atoms = fix_edges(crack_slab)
    qmmm_pot = ForceMixingCarvingCalculator(crack_slab, qm_region_mask,
                                            mm_pot, qm_pot,
                                            buffer_width=buff,
                                            pbc_type=[False, False, True])
    crack_slab.set_calculator(qmmm_pot)

    if args.precon:
        opt = PreconFIRE(crack_slab)
    else:
        opt = FIRE(crack_slab)

    crack_slab.new_array('qm_atoms',qm_region_mask)
    crack_slab.new_array('qm_buffer_atoms',qm_buffer_mask)
    def write_slab(a=crack_slab):
        write_xyz('crack_slab.xyz', a, append=True)
    opt.attach(write_slab)
    opt.run(fmax=args.fmax)
    crack_slab.write('relaxed.xyz')

