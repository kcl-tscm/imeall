import os
import sys

from ase.constraints import FixAtoms
from ase.optimize    import BFGS, FIRE, LBFGS, MDMin, QuasiNewton

from quippy import Atoms, Potential, calc_nye_tensor
from qlab import set_fortran_indexing

from imeall.run_dyn import GBRelax
from imeall.lotf.screw import Dislocation 

set_fortran_indexing=False
gbr = GBRelax()

x = [-1,-1,-1]
y = [ 1,-1, 0]
z = [-1,-1, 2]

s_system = list(x)
s_system.extend(y)
name = ''.join(map(str,s_system)).replace('-','')
name = 'e{0}.xyz'.format(name)

print x,y,z

disloc  = Dislocation(x=x, y=y, z=z, name=name)
disloc.gen_edge_dislocation()
POT_DIR = '/Users/lambert/pymodules/imeall/imeall/potentials' 
eam_pot = os.path.join(POT_DIR,'Fe_Mendelev.xml')

print 'Open Atoms'
at  = Atoms('e{0}.xyz'.format(name))
at.set_cutoff(3.0)
at.calc_connect()


print 'Delete atoms'
at                    = gbr.delete_atoms(grain=at, rcut=1.7)
at.info['OrigHeight'] = at.positions[:,1].max()-at.positions[:,1].min()
r_scale = 1.00894848312
rem = []
for atom in at:
  if atom.position[0] <= 0.716:
    rem.append(atom.index+1)
print 'Removing ', len(rem), ' atoms.'
if len(rem) > 0:
  at.remove_atoms(rem)
else:
  print 'No atoms displaced from unitcell'
pot     = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)
at.set_calculator(pot)
at.write('unrelaxed.xyz')
print 'Running Relaxation'
opt   = FIRE(at)
opt.run(fmax=0.1)

print 'Calculating connectivity and Nye Tensor'

ref_slab  = Atoms('ref_slab.xyz')
ref_slab.set_cutoff(3.0)
ref_slab.calc_connect()

at.set_cutoff(3.0)
at.calc_connect()
alpha     = calc_nye_tensor(at, ref_slab,3,3,at.n)
at.add_property('screw', alpha[2,2,:])
at.add_property('edgex', alpha[2,0,:])
at.add_property('edgey', alpha[2,1,:])
at.write('{0}'.format(name))

