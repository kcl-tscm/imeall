from ase.lattice.cubic import BodyCenteredCubic
from ase.optimize      import BFGS, FIRE, LBFGS, MDMin, QuasiNewton
from qlab              import view
from quippy.io         import AtomsWriter, AtomsReader, write
from quippy.structures import disloc_noam
from quippy import Potential, Atoms, write, supercell, calc_nye_tensor
import numpy as np

#Create unit cell with the orientation:
x = [1,1,-2] 
y = [-1,1,0] 
z = [1,1,1]

screw_slab_unit = BodyCenteredCubic(directions = [x, y, z], size = (1,1,1), symbol='Fe', pbc=(1,1,1),latticeconstant = 2.83)
screw_slab_unit.write('screw_111.xyz')
screw_slab_unit = Atoms(screw_slab_unit)
# set the supercell:
l_x = 150.
l_y = 150.
l_z = 1.

n_x = int(150./screw_slab_unit.get_cell()[0,0])
n_y = int(150./screw_slab_unit.get_cell()[1,1])

print n_x, n_y 
screw_slab_super = supercell(screw_slab_unit, n_x, n_y,1)
#Burgers vector modulus:
b = screw_slab_unit.get_cell()[2,2]
brg_vec = b*np.array([0,0,1])
screw_slab_super.lattice[1,1] += 30.
screw_slab_super.lattice[2,2] += 30.

l = np.array([0,0,1])
screw_slab_super.set_lattice(screw_slab_super.get_cell(), scale_positions=False)
screw_slab_super.write('slab_super.xyz')

super_slab = screw_slab_super.copy()

L_x = screw_slab_super.lattice[1,1]
L_y = screw_slab_super.lattice[2,2]
L_z = screw_slab_super.lattice[3,3]
core = screw_slab_super.pos[:, screw_slab_super.n/2] 

print 'core', core
print 'Dislocation line', l
print 'Number of Atoms in slab', screw_slab_super.n
print 'Burgers Vector', brg_vec

screw_slab_super.set_cutoff(3.0)
screw_slab_super.calc_connect()

disloc_noam(screw_slab_super, core, l, brg_vec) 

screw_slab_super.write('screw_disloc.xyz')
print screw_slab_super.n, super_slab.n
print screw_slab_super.lattice, super_slab.lattice

ref_slab = screw_slab_unit*(n_x,n_y,1)
ref_slab.set_cutoff(3.0)
ref_slab.calc_connect()

alpha = calc_nye_tensor(screw_slab_super, ref_slab, 3,3, super_slab.n)

screw_slab_super.add_property('screw', alpha[3,3,:])
screw_slab_super.add_property('edge', alpha[3,1,:])
screw_slab_super.write('screw_disloc_after.xyz')

#Set Potential
#pot = Potential('IP EAM_ErcolAd ', param_filename='iron_mish.xml')
#pot = Potential('IP EAM_ErcolAd ', param_filename='Fe_Mendelev.xml')
pot = Potential('IP EAM_ErcolAd ', param_filename='Fe_Dudarev.xml')
screw_slab_super.set_calculator(pot)
opt = FIRE(screw_slab_super)
out = AtomsWriter('{0}'.format('{0}_mish_traj.xyz'.format('screw')))

#for i in range(100):
#  opt.run(fmax=0.005, steps=1)
#  out.write(screw_slab_super)
#  if max(np.sum(screw_slab_super.get_forces()**2, axis=1)**0.5) < 0.005:
#    break
#out.close()

opt.run(fmax=0.005)
out.write(screw_slab_super)
ref_slab.set_calculator(pot)

e_bulk = ref_slab.get_potential_energy()/float(len(ref_slab))
e_screw = screw_slab_super.get_potential_energy()/float(len(screw_slab_super))

print 'Energy Screw Dislocation', e_screw
print 'Energy bulk reference', e_bulk
print 'Difference ',e_screw-e_bulk


