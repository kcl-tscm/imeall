import os
import argparse
import numpy as np

from ase.optimize import BFGS, FIRE
from ase.constraints import UnitCellFilter
from ase.lattice.cubic import BodyCenteredCubic
from ase import Atoms as aseAtoms
from fracture.hydrify_cracktips import Hydrify

from quippy import Atoms, Potential, AtomsReader

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--defect', action='store_true')
args = parser.parse_args()

#Hydrogen interatomic distance 0.73
#Hydrogen formation energy -4.
POT_DIR = os.environ['POTDIR']
eam_pot = os.path.join(POT_DIR, 'PotBH.xml')
r_scale = 1.00894848312
pot = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)
alat = 2.83

tetra_pos = alat*np.array([0.5, 0.0, 0.75])
#Structures
gb = BodyCenteredCubic(directions = [[1,0,0], [0,1,0], [0,0,1]],
                       size = (3,3,3), symbol='Fe', pbc=(1,1,1),
                       latticeconstant = alat)

mid_point = 0.5*np.diag(gb.get_cell())
print mid_point
tetra_pos += mid_point
print tetra_pos

gb = Atoms(gb)
if args.defect:
  gb.add_atoms(tetra_pos, 1)
  gb.write('fe_bcc_h.xyz')
else:
  gb.write('fe_bcc.xyz')
