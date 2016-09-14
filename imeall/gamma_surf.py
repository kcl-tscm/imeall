import os
import numpy as np

from   qlab import set_fortran_indexing

import ase.units as units
from   ase.lattice.cubic import Diamond, BodyCenteredCubic
from   ase.optimize      import BFGS, FIRE, LBFGS, MDMin, QuasiNewton
from   ase.constraints   import UnitCellFilter, StrainFilter, FixedLine

from quippy import Atoms, Potential

set_fortran_indexing(False)


#Unit Cell Appropriate
#crack_geom = {'cleavage_plane'    : (1,-1,0),
#              'crack_front'       : (1,1,0),
#              'crack_direction'   : (0,0,1)}
crack_geom = {'cleavage_plane'    : (1,1,0),
              'crack_front'       : (1,-1,0),
              'crack_direction'   : (0,0,1)}

crack_direction = crack_geom['crack_direction']
crack_front     = crack_geom['crack_front']
cleavage_plane  = crack_geom['cleavage_plane']


#fe_unit = BodyCenteredCubic(directions=[crack_direction, cleavage_plane, crack_front],
#                            size=(1,1,1), symbol='Fe', pbc=(1,1,1),
#                            latticeconstant=2.83)
fe_unit = BodyCenteredCubic(directions=[crack_direction, crack_front, cleavage_plane],
                            size=(1,1,1), symbol='Fe', pbc=(1,1,1),
                            latticeconstant=2.83)

#fe_bulk = BodyCenteredCubic(directions=[crack_direction, cleavage_plane, crack_front],
#                            size=(1,8,1), symbol='Fe', pbc=(1,1,1),
#                            latticeconstant=2.83)
nunits  = 8
fe_bulk = BodyCenteredCubic(directions=[crack_direction, crack_front, cleavage_plane],
                            size=(1,1,nunits), symbol='Fe', pbc=(1,1,1),
                            latticeconstant=2.83)

fe_bulk = Atoms(fe_bulk)
fe_unit = Atoms(fe_unit)

ycut    = 5.0 + float(nunits)/2.0*fe_unit.lattice[2,2]
#fe_bulk.center(vacuum=5.0, axis=1)
fe_bulk.center(vacuum=5.0, axis=2)
print 'lattice', fe_bulk.lattice[1,1]
print 'ycut', ycut

#a1 =  fe_bulk.lattice[0,0]
#a2 =  fe_bulk.lattice[2,2]
a1 =  fe_bulk.lattice[0,0]
a2 =  fe_bulk.lattice[1,1]

fe_bulk.write('g001110b0.0.xyz')
POT_DIR = '/Users/lambert/pymodules/imeall/imeall/potentials' 
eam_pot = os.path.join(POT_DIR, 'Fe_Mendelev.xml')
r_scale = 1.00894848312
#eam_pot = os.path.join(POT_DIR, 'iron_mish.xml')
#r_scale = 1.0129007626

pot     = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)
fe_bulk.set_calculator(pot)
print 'Bulk Energy', fe_bulk.get_potential_energy()
refen = fe_bulk.get_potential_energy()
A = fe_bulk.get_volume()/fe_bulk.lattice[2,2]

fil1 = open('ur111110z{0}.dat'.format(nunits), 'w')
fil2 = open('r111110z{0}.dat'.format(nunits), 'w')

for inum, ashift in enumerate(np.arange(0, 1.10, 0.10)):
  fe_shift = fe_bulk.copy()
  fe_shift.set_calculator(pot)
  for at in fe_shift:
    #if at.position[1] > ycut:
    if at.position[2] > ycut:
    #[00-1](110)
      #at.position += ashift*np.array([-a1,0,0])
    #[-1-11](110)
    #  at.position += 0.5*ashift*np.array([a1,0,-a2])
    #[1-11](110)
      at.position += 0.5*ashift*np.array([a1, a2,0])
  print >> fil1, ashift, (units.m**2/units.J)*(fe_shift.get_potential_energy()-refen)/A
  #$strain_mask = [0,0,0,0,0,0]
  #ucf         = UnitCellFilter(fe_shift, strain_mask)
  #fe_shift.rattle(0.05)
  line=[]
  for at in fe_shift:
    #line.append( FixedLine(at.index, (0,1,0)) )
    line.append( FixedLine(at.index, (0,0,1)) )
  fe_shift.set_constraint(line)
  opt   = BFGS(fe_shift)
  opt.run(fmax=0.01, steps=1000)
  fe_shift.write('feb{0}.xyz'.format(ashift))
  print >> fil2, ashift, (units.m**2/units.J)*(fe_shift.get_potential_energy()-refen)/A
fil1.close()
fil2.close()
