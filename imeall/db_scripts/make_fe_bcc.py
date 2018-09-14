import os
import argparse
import numpy as np

from ase.optimize import BFGS, FIRE
from ase.constraints import UnitCellFilter
from ase import Atoms as aseAtoms
from ase.lattice.cubic import BodyCenteredCubic
from quippy import Atoms, Potential, AtomsReader

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--tetra', action='store_true')
parser.add_argument('-o', '--octa', action='store_true')
parser.add_argument('-n', '--supercellsize', nargs='+', type=int, default = [3,3,3])
args = parser.parse_args()

#Hydrogen interatomic distance 0.73
#Hydrogen formation energy -4.
POT_DIR = os.environ['POTDIR']
eam_pot = os.path.join(POT_DIR, 'PotBH.xml')
r_scale = 1.00894848312
pot = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)
#alat = 2.82893
#could just use the proper one as well....
alat = 2.85

sup_cell = args.supercellsize
tetra_pos = alat*np.array([0.5, 0.0, 0.75])
octa_pos = alat*np.array([0.5, 0.5, 0.0])
#Structures
gb = BodyCenteredCubic(directions = [[1,0,0], [0,1,0], [0,0,1]],
                       size = (sup_cell[0],sup_cell[1],sup_cell[2]), symbol='Fe', pbc=(1,1,1),
                       latticeconstant = alat)
mid_point = 0.5*(np.diag(gb.get_cell()))
mid_point = [((sup_cell[0]-1)/2.)*alat for sp in sup_cell]
gb = Atoms(gb)
if args.tetra:
    print 'Tetrahedral Defect'
    tetra_pos += mid_point
    gb.add_atoms(tetra_pos, 1)
    gb.write('bcc_h.xyz')
elif args.octa:
    print 'Octahedral Defect'
    octa_pos += mid_point
    gb.add_atoms(octa_pos, 1)
    gb.write('bcc_h.xyz')
else:
    gb.write('bcc.xyz')
