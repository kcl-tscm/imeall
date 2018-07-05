import re
import os
import sys

import argparse

from matscipy.socketcalc import VaspClient, SocketCalculator
from distutils import spawn

from ase import io
from ase.optimize import FIRE
from ase.optimize.precon import PreconLBFGS
from ase.calculators import vasp
from ase.constraints import FixedLine

from quippy import AtomsWriter, Atoms

from ase.io.xyz import write_xyz

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

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input_file', default=None)
parser.add_argument('-k','--kpar', type=int, default=4)
args = parser.parse_args()

if os.path.isfile('relaxation.traj'):
    print 'using relaxation.traj'
    gam_cell = io.read('relaxation.traj', index='-1')
else:
    print 'using init_relaxed.xyz'
    gam_cell = io.read('init_relaxed.xyz')

magmoms = [len(gam_cell), 2.24]

vasp_args = dict(xc='PBE', amix=0.22, amin=0.02, bmix=0.6, amix_mag=1.0, bmix_mag=0.9, lorbit=11,
                 kpts=[1, 6, 6], kpar=args.kpar, lreal='auto', nelmdl=-15, ispin=2, prec='Accurate', ediff=1.e-4,
                 encut=480, nelm=100, algo='VeryFast', lplane=False, lwave=False, lcharg=False, istart=0,
                 magmom=magmoms, maxmix=25, #https://www.vasp.at/vasp-workshop/slides/handsonIV.pdf #for badly behaved clusters.
                 voskown=0, ismear=1, sigma=0.1, isym=2, iwavpr=11) 

mpirun = spawn.find_executable('mpirun')
vasp = '/home/mmm0007/vasp/vasp.5.4.1/bin/vasp_std'

vasp_client = VaspClient(client_id=0, npj=96, ppn=1,
                         exe=vasp, mpirun=mpirun, parmode='mpi',
                         ibrion=13, nsw=1000000,
                         npar=6, **vasp_args)

traj = io.Trajectory('relaxation.traj','a', gam_cell)

qm_pot = SocketCalculator(vasp_client)
gam_cell.set_calculator(qm_pot)

fixed_line=[]
for at in gam_cell:
    fixed_line.append(FixedLine(at.index, (1,0,0)))
gam_cell.set_constraint(fixed_line)
opt = PreconLBFGS(Atoms(gam_cell))
#opt.attach(traj_writer, interval=1)
opt.attach(traj.write, interval=1)
opt.run(fmax=0.01)
traj.close()

