import os
import numpy as np
from ase.lattice.cubic import BodyCenteredCubic
from ase.optimize      import BFGS, FIRE, LBFGS, MDMin, QuasiNewton

from qlab              import view
from quippy.io         import AtomsWriter, AtomsReader, write
from quippy.structures import disloc_noam
from quippy            import Potential, Atoms, write, supercell, calc_nye_tensor
from quippy            import set_fortran_indexing
from quippy.farray     import fzeros, frange, unravel_index, farray

#Create unit cell with the orientation:

class Dislocation(object):
  def __init__(self, x=[1,1,-2], y=[-1,1,0],z=[1,1,1], name='defect.xyz'):
    self.x = x 
    self.y = y 
    self.z = z
    self.l_x =  20.
    self.l_y =  20.
    self.l_z =  20.
    self.name = name

  def gen_screw_dislocation(self):
    screw_slab_unit = BodyCenteredCubic(directions = [self.x, self.y, self.z], 
                                        size = (1,1,1), symbol='Fe', pbc=(1,1,1),
                                        latticeconstant = 2.83)
    screw_slab_unit = Atoms(screw_slab_unit)

# set the supercell:
    n_x = int(self.l_x/screw_slab_unit.get_cell()[0,0])
    n_y = int(self.l_y/screw_slab_unit.get_cell()[1,1])
    n_z = 2

    screw_slab_super = supercell(screw_slab_unit, n_x, n_y, n_z)
    ref_slab = screw_slab_unit*(n_x,n_y,n_z)
    ref_slab.write('ref_slab.xyz')
#Burgers vector modulus:
    b = screw_slab_unit.get_cell()[2,2]
    brg_vec = b*np.array([0,0,1])
    print b

    vacuum = 70.
    screw_slab_super.lattice[0,0] += vacuum
    screw_slab_super.lattice[1,1] += vacuum

    disloc_l   = np.array([0,0,1])
    screw_slab_super.set_lattice(screw_slab_super.get_cell(), scale_positions=False)
    super_slab = screw_slab_super.copy()

    L_x  = screw_slab_super.lattice[0,0]
    L_y  = screw_slab_super.lattice[1,1]
    L_z  = screw_slab_super.lattice[2,2]
    core =  np.array([(L_x-vacuum)/2., (L_y-vacuum)/2., L_z/2.])

    screw_slab_super.set_cutoff(3.0)
    screw_slab_super.calc_connect()
    disloc_noam(screw_slab_super, core, disloc_l, brg_vec) 
    screw_slab_super.info['core'] = core 
    screw_slab_super.write('s{0}.xyz'.format(self.name))


  def gen_edge_dislocation(self):
    screw_slab_unit = BodyCenteredCubic(directions = [self.x, self.y, self.z], 
                                        size = (1,1,1), symbol='Fe', pbc=(1,1,1),
                                        latticeconstant = 2.83)
    screw_slab_unit = Atoms(screw_slab_unit)

    n_x = int(self.l_x/screw_slab_unit.get_cell()[0,0])
    n_y = int(self.l_y/screw_slab_unit.get_cell()[1,1])
    n_z = 1

    screw_slab_super = supercell(screw_slab_unit, n_x, n_y, n_z)
    screw_slab_unit.write('ref_slab.xyz')

    b = screw_slab_unit.get_cell()[0,0]
    brg_vec = b*np.array([1,0,0])

    vacuum = 70.
    screw_slab_super.lattice[1,1] += vacuum

    disloc_l = np.array([0,0,1])
    screw_slab_super.set_lattice(screw_slab_super.get_cell(), scale_positions=False)
    super_slab = screw_slab_super.copy()

    L_x  = screw_slab_super.lattice[0,0]
    L_y  = screw_slab_super.lattice[1,1]
    L_z  = screw_slab_super.lattice[2,2]
#Maintain Periodicity along x and z for an edge dislocation:
    core =  np.array([(L_x)/2., (L_y-vacuum)/2., (L_z)/2.])

    screw_slab_super.set_cutoff(3.0)
    screw_slab_super.calc_connect()
    disloc_noam(screw_slab_super, core, disloc_l, brg_vec) 
    screw_slab_super.info['core'] = core 
    screw_slab_super.write('e{0}.xyz'.format(self.name))
    
