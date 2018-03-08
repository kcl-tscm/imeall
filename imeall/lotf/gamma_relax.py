import re
import sys

import os
from matscipy.socketcalc import VaspClient, SocketCalculator
from distutils import spawn

from ase import io
from ase.optimize import FIRE
from ase.optimize.precon import PreconLBFGS
from ase.calculators import vasp
from ase.constraints import FixedLine

from quippy import AtomsWriter, Atoms

def read_vasp_args():
    kv_pairs = re.compile(r'(.*?)=(.*?)\n')
    vasp_args = {}
    with open('INCAR', 'r') as f:
        vasp_incar = f.read()
    kv_pairs = kv_pairs.findall(vasp_incar)
    print kv_pairs
    for kv in kv_pairs:
        vasp_args[kv[0].strip().lower()] = kv[1].strip().lower()
    return vasp_args

#vasp_pot = vasp.Vasp(restart=False)
#vasp_pot.read_incar()
#vasp_pot.read_kpoints()
#vasp_pot.read_potcar()
#print vasp_pot.input_params
#vasp_args = read_vasp_args()
#for k,v in vasp_args.items():
#    print k, v
#vasp_pot.set(**vasp_args)

gam_cell = io.read('init_relaxed.xyz')
magmoms = [len(gam_cell), 2.24]

vasp_args = dict(xc='PBE', amix=0.22, amin=0.02, bmix=0.9, amix_mag=1.1, bmix_mag=1.0,
                 kpts=[3, 3, 1], kpar=8, lreal='auto', nelmdl=-15, ispin=2, prec='Accurate', ediffg=-1.e-4,
                 encut=420, nelm=120, algo='VeryFast', lplane=False, lwave=False, lcharg=False, istart=0,
                 magmom=magmoms, maxmix=-25, #https://www.vasp.at/vasp-workshop/slides/handsonIV.pdf #for badly behaved clusters.
                 voskown=0, ismear=1, sigma=0.1, isym=2) # possibly try iwavpr=12, should be faster if it works


mpirun = spawn.find_executable('mpirun')
vasp = '/home/mmm0007/vasp/vasp.5.4.1/bin/vasp_std'

vasp_client = VaspClient(client_id=0, npj=48, ppn=1,
                         exe=vasp, mpirun=mpirun, parmode='mpi',
                         ibrion=13, nsw=1000000,
                         npar=6, **vasp_args)

def pass_trajectory_context(trajectory, dynamics):
    def traj_writer(dynamics):
        trajectory.write(dynamics.atoms)
    return traj_writer

trajectory = AtomsWriter('relaxation.xyz')

qm_pot = SocketCalculator(vasp_client)
gam_cell.set_calculator(qm_pot)

#fixed_line=[]
#for at in gam_cell:
#    fixed_line.append(FixedLine(at.index, (1,0,0)))
#gam_cell.set_constraint(fixed_line)

opt = PreconLBFGS(Atoms(gam_cell))
opt.attach(pass_trajectory_context(trajectory, opt), 1, opt)
opt.run(fmax=0.01)
