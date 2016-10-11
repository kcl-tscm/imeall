from quippy import bcc
from quippy.potential import Potential, Minim
import numpy as np

fe_bulk = bcc(2.83)
fe_bulk.set_atoms(26)
#center_cell = np.diag(fe_bulk.get_cell())/4.
#fe_bulk.add_atoms(center_cell,1)
#fe_bulk.set_atoms(fe_bulk.Z)

eam_pot = 'iron_mish.xml'
eam_pot = 'Fe_Mendelev.xml'
eam_pot = 'quip.xml'
eam_pot = 'Fe_Ackland.xml'
eam_pot = 'PotBH.xml'
print eam_pot

r_scale = 1.0
pot     = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale=1.0', r_scale=r_scale, param_filename=eam_pot)
fe_bulk.set_calculator(pot)
print pot.cutoff()
initial_a0 = fe_bulk.get_cell()[0][0]

print 'Before Relaxation'
print 'Energy', fe_bulk.get_potential_energy()/len(fe_bulk)
print fe_bulk.get_cell().round(5)

minim   = Minim(fe_bulk, relax_positions=True, relax_cell=True)
minim.run(fmax=1e-4)

print 'After Relaxation'
print 'Energy', fe_bulk.get_potential_energy()/len(fe_bulk)
print fe_bulk.get_cell().round(5)

fe_bulk.write('FeH.xyz')

final_a0 = fe_bulk.get_cell()[0][0]

print 'Rescaling Potential: '
r_scale = final_a0/initial_a0
print r_scale

fe_bulk = bcc(2.83)
fe_bulk.set_atoms(26)
print 'IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale)
pot     = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)
fe_bulk.set_calculator(pot)
print fe_bulk.get_potential_energy()/len(fe_bulk)

minim   = Minim(fe_bulk, relax_positions=True, relax_cell=True)
minim.run(fmax=1e-4)

print 'After rescaled Relaxation'
print 'Energy', fe_bulk.get_potential_energy()/len(fe_bulk)
print fe_bulk.get_cell().round(5)

fe_bulk.write('FeH.xyz')


