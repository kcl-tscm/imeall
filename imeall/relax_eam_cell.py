from quippy.potential import Potential, Minim
from quippy import bcc

fe_bulk = bcc(2.83)
fe_bulk.set_atoms(26)

eam_pot = 'iron_mish.xml'
eam_pot  = 'Fe_Ackland.xml'
eam_pot = 'Fe_Mendelev.xml'
print eam_pot

r_scale = 1.0
pot     = Potential('IP EAM_ErcolAd do_rescale=T', r_scale=r_scale, param_filename=eam_pot)
fe_bulk.set_calculator(pot)
print pot.cutoff()
initial_a0 = fe_bulk.get_cell()[0][0]

print 'Before Relaxation'
print fe_bulk.get_potential_energy()/len(fe_bulk)
print fe_bulk.get_cell().round(2)

minim   = Minim(fe_bulk, relax_positions=True, relax_cell=True)
minim.run(fmax=1e-4)

print 'After Relaxation'
print fe_bulk.get_potential_energy()/len(fe_bulk)
print fe_bulk.get_cell().round(2)

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
print fe_bulk.get_potential_energy()/len(fe_bulk)
print fe_bulk.get_cell().round(2)

