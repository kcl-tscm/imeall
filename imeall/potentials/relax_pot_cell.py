import sys
import argparse
import numpy as np
from quippy import bcc
from quippy.potential import Potential, Minim



parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file', help='name of input file')
parser.add_argument('-pt', '--pot_type', help='Potential type, GAP. EAM.')
args = parser.parse_args()

fe_bulk = bcc(2.83)
fe_bulk.set_atoms(26)

print 'Input File', args.input_file

r_scale = 1.0
if args.pot_type =='EAM':
  pot = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale=1.0', r_scale=r_scale, param_filename=args.input_file)
elif args.pot_type =='GAP':
  pot = Potential('IP GAP', param_filename = args.input_file)
else:
  sys.exit('Invalid pot type.')

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

print 'IP {1} do_rescale_r=T r_scale={0}'.format(r_scale, args.pot_type)
if args.pot_type == 'EAM':
  pot     = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=args.input_file)
elif args.pot_type == 'GAP':
  #rescale not yet implemented
  #pot     = Potential('IP GAP do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=args.input_file)
  pot     = Potential('IP GAP'.format(r_scale), param_filename=args.input_file)
else:
  sys.exit('Invalid pot_type')

fe_bulk.set_calculator(pot)
print fe_bulk.get_potential_energy()/len(fe_bulk)

minim   = Minim(fe_bulk, relax_positions=True, relax_cell=True)
minim.run(fmax=1e-4)

print 'After rescaled Relaxation'
print 'Energy', fe_bulk.get_potential_energy()/len(fe_bulk)
print fe_bulk.get_cell().round(5)
fe_bulk.write('Fe_rescaled_relax.xyz')


