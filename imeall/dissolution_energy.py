import os
import numpy as np

from ase.optimize import BFGS, FIRE
from ase.constraints import UnitCellFilter
from ase.lattice.cubic import BodyCenteredCubic
from ase import Atoms as aseAtoms
from fracture.hydrify_cracktips import Hydrify

from quippy import Atoms, Potential, AtomsReader

POT_DIR = os.environ['POTDIR']
eam_pot = os.path.join(POT_DIR, 'PotBH.xml')
r_scale = 1.00894848312
pot = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)
alat = 2.83
#Structures
gb = BodyCenteredCubic(directions = [[1,0,0], [0,1,0], [0,0,1]],
                       size = (4,4,4), symbol='Fe', pbc=(1,1,1),
                       latticeconstant = alat)


tetra_pos = alat*np.array([0.25, 0.0, 0.5])

h2 = aseAtoms('H2', positions=[[0, 0, 0],[0, 0, 0.7]])
h2 = Atoms(h2)

gb = Atoms(gb)
gb_h = gb.copy()
gb_h.add_atoms(tetra_pos, 1)


#caclulators
gb.set_calculator(pot)
h2.set_calculator(pot)
gb_h.set_calculator(pot)
gb_h.write('hydrogen_bcc.xyz')

#Calc Hydrogen molecule energy
opt = BFGS(h2)
opt.run(fmax=0.0001)
E_h2  = h2.get_potential_energy()


strain_mask  = [1,1,1,0,0,0]
ucf = UnitCellFilter(gb_h, strain_mask)

#opt = BFGS(gb_h)
opt = FIRE(ucf)
opt.run(fmax=0.0001)
E_gb  = gb.get_potential_energy()
E_gbh = gb_h.get_potential_energy()
E_dis = E_gbh - E_gb - 0.5*E_h2

print 'E_gb', E_gb
print 'E_gbh', E_gbh
print 'H2 Formation Energy', E_h2
print 'Dissolution Energy', E_dis

