import os
import numpy as np
from ase.lattice.cubic   import BodyCenteredCubic
from quippy       import Atoms, Potential, AtomsReader
from ase.optimize import BFGS
from fracture.hydrify_cracktips import Hydrify
from ase import Atoms as aseAtoms

POT_DIR = os.environ['POTDIR']
eam_pot = os.path.join(POT_DIR, 'PotBH.xml')
r_scale = 1.00894848312
pot     = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)

#Structures
gb      = BodyCenteredCubic(directions = [[1,0,0], [0,1,0], [0,0,1]],
                            size = (8,8,8), symbol='Fe', pbc=(1,1,1),
                            latticeconstant = 2.83)

calc_elast_dipole=True
if calc_elast_dipole:
  cell = gb.get_cell()
  cell = cell*1.01
  gb.set_cell(cell, scale_atoms=True)

print cell[0,0], cell[1,1], cell[2,2]
print cell[0,0]*0.01, cell[1,1]*0.01, cell[2,2]*0.01
print cell[0,0]*cell[1,1]*cell[2,2]
1/0

tetra_pos = 2.83*np.array([0.25, 0.0, 0.5])

h2 = aseAtoms('H2', positions=[[0, 0, 0],[0, 0, 0.7]])
h2 = Atoms(h2)
gb = Atoms(gb)
gb_h = gb.copy()
gb_h.add_atoms(tetra_pos, 1)

#caclulators
gb.set_calculator(pot)
h2.set_calculator(pot)
gb_h.set_calculator(pot)

opt = BFGS(h2)
opt.run(fmax=0.001)
E_h2  = h2.get_potential_energy()

opt = BFGS(gb_h)
opt.run(fmax=0.001)
E_gb  = gb.get_potential_energy()
E_gbh = gb_h.get_potential_energy()
E_dis = E_gbh - E_gb - 0.5*E_h2

print 'E_gb', E_gb
print 'E_gbh', E_gbh
print 'H2 Formation Energy', E_h2
print 'Dissolution Energy', E_dis



