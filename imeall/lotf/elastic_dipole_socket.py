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

def gen_impurity(symbol='H', tetrahedral=True, sup_cell= [5,5,5]):
    #Molybdenum is 42!!!
    alat = 2.82893
    atomic_number_dict = {'H':1, 'B':5, 'C':6, 'P':15, 'Mn':25, 'Mo':42}

    tetra_pos = alat*np.array([0.5, 0.0, 0.75])
    octa_pos = alat*np.array([0.5, 0.5, 0.0])
    #Structures
    ats = BodyCenteredCubic(directions = [[1,0,0], [0,1,0], [0,0,1]],
                            size = (sup_cell[0],sup_cell[1],sup_cell[2]), symbol='Fe', pbc=(1,1,1),
                            latticeconstant = alat)

    mid_point = 0.5*(np.diag(gb.get_cell()))
    mid_point = [((sup_cell[0]-1)/2.)*alat for sp in sup_cell]

    if tetrahedral:
        tetra_pos += mid_point
        ats.append(Atom(position=tetra_pos, symbol=symbol))
    else:
        octa_pos += mid_point
        ats.append(Atom(position=octa_pos, symbol=symbol))

    return ats


parser = argparse.ArgumentParser()
parser.add_argument('-i','--input_file', default=None)
parser.add_argument('-k','--kpar', type=int, default=4)
parser.add_argument('-s','--symbol', default='H', help="Chemical symbol for impurity introduced")
parser.add_argument('-o','--octahedral', help="Add impurity at tetrahedral position in cell.", action='store_true')
args = parser.parse_args()

gam_cell = gen_impurity(symbol = args.symbol, tetrahedral=not args.octahedral)
gam_cell.write('init_defect_cell.xyz')

magmoms = [len(gam_cell), 2.24]
vasp_args = dict(xc='PBE', amix=0.22, amin=0.02, bmix=0.6, amix_mag=1.0, bmix_mag=0.9,
                 kpts=[4, 4, 4], kpar=args.kpar, lreal='auto', nelmdl=-15, ispin=2, prec='Accurate', ediff=1.e-5,
                 encut=480, nelm=100, algo='VeryFast', lplane=False, lwave=False, lcharg=False, istart=0,
                 magmom=magmoms, maxmix=25, #https://www.vasp.at/vasp-workshop/slides/handsonIV.pdf #for badly behaved clusters.
                 voskown=0, ismear=1, sigma=0.1, isym=2, iwavpr=11) 

mpirun = spawn.find_executable('mpirun')
vasp = '/home/mmm0007/vasp/vasp.5.4.1/bin/vasp_std'

vasp_client = VaspClient(client_id=0, npj=96, ppn=1,
                         exe=vasp, mpirun=mpirun, parmode='mpi',
                         ibrion=13, nsw=1000000,
                         npar=6, **vasp_args)

traj = io.Trajectory('relaxation.traj','w', gam_cell)

qm_pot = SocketCalculator(vasp_client)
gam_cell.set_calculator(qm_pot)

fixed_line=[]
opt = PreconLBFGS(Atoms(gam_cell))
opt.attach(traj.write, interval=1)
opt.run(fmax=0.01)
traj.close()

#remove defect atom with symbol 
del gam_cell[[atom.symbol == args.symbol for atom in gam_cell]]
#compute induced dipole forces.
gam_cell.get_forces()
gam_cell.write('forces.xyz')

