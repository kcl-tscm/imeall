import re
import os
import sys
import argparse
import numpy as np

from ase import io, Atoms, Atom
from ase.optimize import FIRE
from ase.calculators import vasp
from ase.io.xyz import write_xyz
from ase.constraints import FixedLine
from ase.optimize.precon import PreconLBFGS
from ase.lattice.cubic import BodyCenteredCubic

from distutils import spawn

from matscipy.socketcalc import VaspClient, SocketCalculator


def gen_segrant_dirs(template_dir="../templates", template_name="ed_calc.sh", tetrahedral=False):
    """
    Create directories for all segregant/elastic dipole calculations. And copy
    submission templates. (This routine can me modified for different segregants 
    and requires a submission script in a templates directory folder).

    Args:
        template_dir(str): location of directory containing template submission script.
        template_name(str): name of template submission script.
    """
    segregants = ['H', 'C', 'N', 'P', 'Mn', 'Mo', 'B']
    for num, segregant in enumerate(segregants):
        target_dir = "{}_Tet".format(segregant)
        try:
            os.mkdir(target_dir)
        except OSError:
            print 'Dir already exists.'
        with pushd(target_dir) as ctx0:
            template_pbs = os.path.join(template_dir, template_name)
            shutil.copy(template_pbs, './')
            with open('ed_calc.sh','r') as f:
                pbs_file =f.read()
            if tetrahedral: 
                tet_flag = '-t'
            else:
                tet_flag = ''
            pbs_file = pbs_file.format(jobname=target_dir, symbol=segregant, tet_flag=tet_flag)
            with open('ed_calc.sh','w') as f:
                print >> f, pbs_file

def gen_impurity(symbol='H', tetrahedral=True, sup_cell= [5,5,5]):
    """
    Add an element of type symbol to a supercell defined by the size sup_cell.
    If tetrahedral is true the element is initially placed in a tetrahedral position
    otherwise it is placed in an octahedral position.
    """
    #Molybdenum is 42!!!
    alat = 2.82893
    atomic_number_dict = {'H':1, 'B':5, 'C':6, 'P':15, 'Mn':25, 'Mo':42}

    tetra_pos = alat*np.array([0.5, 0.0, 0.75])
    octa_pos = alat*np.array([0.5, 0.5, 0.0])
    #Structures
    ats = BodyCenteredCubic(directions = [[1,0,0], [0,1,0], [0,0,1]],
                            size = (sup_cell[0],sup_cell[1],sup_cell[2]), symbol='Fe', pbc=(1,1,1),
                            latticeconstant = alat)

    mid_point = 0.5*(np.diag(ats.get_cell()))
    mid_point = [((sup_cell[0]-1)/2.)*alat for sp in sup_cell]

    if tetrahedral:
        tetra_pos += mid_point
        ats.append(Atom(position=tetra_pos, symbol=symbol))
    else:
        octa_pos += mid_point
        ats.append(Atom(position=octa_pos, symbol=symbol))
    return ats

parser = argparse.ArgumentParser()
parser.add_argument('-s','--symbol', default='H', help="Chemical symbol for impurity introduced")
parser.add_argument('-o','--octahedral', help="Add impurity at octahedral position in cell.", action='store_true')
parser.add_argument('-k','--kpar', type=int, default=6)
parser.add_argument('-f','--fmax', help='max force relaxation', type=float, default=0.01)
parser.add_argument('-r','--restart', help="restart", action='store_true')
parser.add_argument('-n','--no_relax', help="Do not relax cell, just delete impurity from the relaxation.traj file and calc forces.",  action='store_true')
args = parser.parse_args()

tetrahedral = not args.octahedral
if not args.restart:
    gam_cell = gen_impurity(symbol = args.symbol, tetrahedral=tetrahedral)
    gam_cell.write('init_defect_cell.xyz')
else:
    print 'restarting from relaxation.traj'
    gam_cell = io.read('relaxation.traj', index='-1')


magmoms = [len(gam_cell), 2.24]
vasp_args = dict(xc='PBE', amix=0.22, amin=0.02, bmix=0.6, amix_mag=1.0, bmix_mag=0.9, symprec=1e-8,
                 kpts=[4, 4, 4], kpar=args.kpar, lreal='auto', nelmdl=-15, ispin=2, prec='Accurate', ediff=1.e-5,
                 encut=480, nelm=100, algo='VeryFast', lplane=False, lwave=False, lcharg=False, istart=0,
                 magmom=magmoms, maxmix=35, #https://www.vasp.at/vasp-workshop/slides/handsonIV.pdf #for badly behaved clusters.
                 voskown=0, ismear=1, sigma=0.1, isym=2, iwavpr=11) 

mpirun = spawn.find_executable('mpirun')
vasp = '/home/mmm0007/vasp/vasp.5.4.1/bin/vasp_std'

vasp_client = VaspClient(client_id=0, npj=96, ppn=1,
                         exe=vasp, mpirun=mpirun, parmode='mpi',
                         ibrion=13, nsw=1000000,
                         npar=4, **vasp_args)

if not args.no_relax:
    traj = io.Trajectory('relaxation.traj','a', gam_cell)
    qm_pot = SocketCalculator(vasp_client)
    gam_cell.set_calculator(qm_pot)
    opt = PreconLBFGS(gam_cell)
    opt.attach(traj.write, interval=1)
    opt.run(fmax=args.fmax)
    traj.close()
    qm_pot.shutdown()

#remove defect atom with symbol 
del gam_cell[[atom.symbol == args.symbol for atom in gam_cell]]
defect_cell = gam_cell.copy()
defect_cell.write('no_impurity.xyz')

#Need accurate forces
vasp_args['ediff']=1e-5
vasp_client = VaspClient(client_id=0, npj=96, ppn=1,
                         exe=vasp, mpirun=mpirun, parmode='mpi',
                         ibrion=13, nsw=1000000,
                         npar=4, **vasp_args)

qm_pot = SocketCalculator(vasp_client)
#compute induced dipole forces.
defect_cell.set_calculator(qm_pot)
defect_cell.get_forces()
defect_cell.write('forces.xyz')
qm_pot.shutdown()

