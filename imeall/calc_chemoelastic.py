import os
import sys
import argparse
import subprocess
import numpy as np
import ase.units as units
from quippy import Atoms, Potential, AtomsReader, set_fortran_indexing
from imeall.models import PotentialParameters

set_fortran_indexing(False)
def strain_energy(ats):
  cell = ats.get_cell()
  A = cell[0][0]*cell[1][1]
  z_height = cell[2][2]
  ener_z = np.transpose(np.vstack((ats.get_potential_energies(), [at.position[2] for at in ats])))
  ener_z = np.array(sorted(ener_z, key=lambda x: x[1]))
  cursor = 0.2
  elastic_energy = []
  while cursor < (z_height + 0.2):
    try:
      cum_energy = 16.02*sum([x+4.01298214176 for x in np.array(filter(lambda x: x[1]<= cursor, ener_z))[:,0]])/(A)
    except IndexError:
      cursor += 0.2 #initial step doesn't capture atoms
      continue 
    elastic_energy.append((cursor, cum_energy))
    cursor += 0.2
  return elastic_energy

def convolution(original_curve, z_height=79.0):
  sigma = 0.2
  norm = np.sqrt(2*np.pi)
  dx = float(z_height-0.2)/len(original_curve)
  gx = np.arange(-3*sigma, 3*sigma, dx)
  gaussian = np.exp(-(gx/sigma)**2/2.0)/norm
  result = np.convolve(original_curve, gaussian, mode="full")
  return result

def calc_chemomechanical(ats):
  """
  method:`chemoelastic_percentage` requires atoms object to have at least structure_type
          and local_energy property. bcc lattice structure_type=3.
  """
#case quip types to numpy arrays stack and transpose
  loc_en = np.array(ats.properties['local_energy'])
  struct_type = np.array(ats.properties['structure_type'])
  joined = np.vstack((loc_en, struct_type)).transpose()
  cell = ats.get_cell()
  A = cell[0][0]*cell[1][1]
#compute relative total energy contributions of the
#two types
  gs = np.zeros(np.shape(joined))
#zero energies
  gs[:,0] = -4.01298
  joined = joined-gs
#select bcc and non_bcc
  non_bcc = joined[joined[:,1]==3]
  bcc = joined[joined[:,1] != 3]
  chemical_energy = np.sum(non_bcc[:,0])
  elastic_energy = np.sum(bcc[:,0])
  total_energy = chemical_energy + elastic_energy
  gb_energy = 16.02*(total_energy)/(2.*A)
  print 'Chemical {}%, Elastic {}%'.format(round(100*chemical_energy/total_energy, 2), 
                                           round(100*elastic_energy/total_energy, 2))
  return [(chemical_energy/total_energy)*gb_energy, (elastic_energy/total_energy)*gb_energy, gb_energy]

def calc_chemoelast(input_file):
  potparam = PotentialParameters()
  ener_bulk_dict = potparam.gs_ener_per_atom()
  r_scale_dict = potparam.eam_rscale()
  r_scale = r_scale_dict['PotBH.xml']
  E_bulk = ener_bulk_dict['PotBH.xml']
  try:
    POT_DIR = os.environ ['POTDIR']
  except KeyError:
    sys.exit("PLEASE SET export POTDIR='path/to/potfiles/'")
  eam_pot = 'PotBH.xml'
  eam_pot = os.path.join(POT_DIR, eam_pot)
  pot = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)
  ats = AtomsReader(input_file)[-1]
  ats.set_calculator(pot)
  gb_energy = potparam.calc_e_gb(ats, E_bulk)
  print gb_energy, 'J/^2m'
  ats.write('full.xyz')
  elastic_energy = strain_energy(ats)
  with open('elast.dat','w') as f:
    for x in elastic_energy:
      print >> f, x[0], x[1]
#generates output.xyz
  args_str =  'ovitos /Users/lambert/pymodules/imeall/imeall/ovito_scripts/attach_cna.py -i {input_file}'.format(input_file='full.xyz').split()
  job = subprocess.Popen(args_str)
  job.wait()
  ats = Atoms('output.xyz')
#print the three contributions
  x = calc_chemomechanical(ats)
  try:
    assert round(gb_energy,2) == round(x[2],2)
  except AssertionError:
    print "WARNING ENERGIES DON'T MATCH", gb_energy, x[2]
  return x

if __name__=='__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--input_file', '-i', help='Input grain boundary struct file to produce cumulative energy')
  args = parser.parse_args()
  calc_chemoelast(input_file=args.input_file)
