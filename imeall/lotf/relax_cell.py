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
from argparse import ArgumentParser

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

parser = ArgumentParser()

parser.add_argument('-kp', '--kpar', help='kpoints par', default = 8, type=int)
parser.add_argument('-kx', '--kx', help='kpoints x', default = 1, type=int)
parser.add_argument('-ky', '--ky', help='kpoints y', default = 2, type=int)
parser.add_argument('-kz', '--kz', help='kpoints z', default = 16, type=int)
parser.add_argument('-o', '--outcar', help='read OUTCAR file instead of relaxation', action='store_true')
parser.add_argument('-p', '--poscar', help='read POSCAR file instead of relaxation', action='store_true')
parser.add_argument('-f', '--fmax', help='max force eV/A', default = 0.01, type=float)
parser.add_argument('-s', '--isym', help='isym_flag', default=2, type=int)
parser.add_argument('-n', '--npj', help='num processors', default=48, type=int)
args = parser.parse_args()

if args.outcar:
    print 'READING OUTCAR'
    gam_cell = io.read('OUTCAR', index='-1')
elif args.poscar:
    print 'READING POSCAR'
    gam_cell = io.read('POSCAR', index='-1')
else:
    if os.path.isfile('./relaxation.traj'):
        print 'Starting from relaxation.traj'
        gam_cell = io.read('relaxation.traj', index='-1')
    else:
        print 'Starting from init_relaxed.xyz'
        gam_cell = io.read('init_relaxed.xyz')

magmoms = [len(gam_cell), 2.24]

vasp_args = dict(xc='PBE', amix=0.18, amin=0.018, bmix=0.6, amix_mag=1.0, bmix_mag=0.8,
                 kpts=[args.kx, args.ky, args.kz], kpar=args.kpar, lreal='auto', nelmdl=-15, ispin=2, prec='Accurate', ediff=1.e-4,
                 encut=420, nelm=100, algo='VeryFast', lplane=False, lwave=False, lcharg=False, istart=0, iwavpr=11,
                 magmom=magmoms, maxmix=25, #https://www.vasp.at/vasp-workshop/slides/handsonIV.pdf #for badly behaved clusters.
                 voskown=0, ismear=1, sigma=0.1, isym=args.isym) # possibly try iwavpr=12, should be faster if it works

mpirun = spawn.find_executable('mpirun')
vasp = '/home/mmm0007/vasp/vasp.5.4.1/bin/vasp_std'

vasp_client = VaspClient(client_id=0, npj=args.npj, ppn=1,
                         exe=vasp, mpirun=mpirun, parmode='mpi',
                         ibrion=13, nsw=1000000,
                         npar=6, **vasp_args)

#def pass_trajectory_context(trajectory, dynamics):
#    def traj_writer(dynamics):
#        trajectory.write(dynamics.atoms)
#    return traj_writer
#trajectory = AtomsWriter('relaxation.xyz')

traj = io.Trajectory('relaxation.traj','w', gam_cell)

qm_pot = SocketCalculator(vasp_client)
gam_cell.set_calculator(qm_pot)

opt = PreconLBFGS(gam_cell)
opt.attach(traj.write, interval=1)
opt.run(fmax=args.fmax)

traj.close()
