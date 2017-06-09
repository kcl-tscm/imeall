import os
import sys
import ase.units as units
import numpy as np
from quippy import Atoms, Potential, AtomsReader
from imeall.models import PotentialParameters

def strain_energy(ats):
  cell = ats.get_cell()
  A = cell[0][0]*cell[1][1]
  z_height = cell[2][2]
  ener_z = np.transpose(np.vstack((ats.get_potential_energies(), [at.position[2] for at in ats])))
  ener_z = np.array(sorted(ener_z, key=lambda x: x[1]))

  cursor = 0.2
  elastic_energy = []
  while cursor < (z_height + 0.2):
    cum_energy = 16.02*sum([x+4.01298214176 for x in np.array(filter(lambda x: x[1]<= cursor, ener_z))[:,0]])/(A)
    elastic_energy.append((cursor, cum_energy))
    cursor += 0.2
  return elastic_energy


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

ats = AtomsReader('0011626170_v6bxv2_tv0.3bxv0.4_d1.6z_traj.xyz')[-1]
ats.set_calculator(pot)
print potparam.calc_e_gb(ats, E_bulk), 'J/^2m'
ats.write('full.xyz')

elastic_energy = strain_energy(ats)  
for x in elastic_energy:
  print x[0], x[1]

