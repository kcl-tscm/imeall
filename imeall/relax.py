import os
import sys
import json
import shutil
import ase.io        
from   cStringIO           import StringIO
from   ase.optimize.sciopt import SciPyFminCG
from   quippy              import Atoms, Potential, frange
from   ase.constraints     import UnitCellFilter, StrainFilter
from   quippy.io           import AtomsWriter, AtomsReader, write
from   ase.optimize        import BFGS, FIRE, LBFGS, MDMin, QuasiNewton
import numpy as np
from pprint import pprint

def converged(grain, smax, fmax):
  maxstress = max(grain.get_stress().ravel())
  rmsforces = np.sum(grain.get_forces()**2, axis=1)**0.5
  maxforce  = max(rmsforces)
  if maxforce < fmax and maxstress < smax:
    return True
  return False


with open('subgb.json', 'r') as outfile:
  j_dict = json.load(outfile)

try: 
  param_file = j_dict['param_file']
  if param_file == 'iron_mish.xml':
    eam_pot = '/users/k1511981/sharedscratch/grain_boundaries/iron_mish.xml'
  elif param_file == 'Fe_Mendelev.xml':
    eam_pot = '/users/k1511981/sharedscratch/grain_boundaries/Fe_Mendelev.xml'
except KeyError:
  print 'No EAM pot specified defaulting to Mendelev'
  eam_pot = '/users/k1511981/sharedscratch/grain_boundaries/Fe_Mendelev.xml'

eam_pot = '/users/k1511981/sharedscratch/grain_boundaries/Fe_Mendelev.xml'

print eam_pot

pot_file = eam_pot.split('/')[-1]
grain   = Atoms('{0}'.format(sys.argv[1]))
pot     = Potential('IP EAM_ErcolAd', param_filename=eam_pot)
grain.set_calculator(pot)
E_gb_init   = grain.get_potential_energy()

alpha       = E_gb_init
out         = AtomsWriter('{0}'.format('{0}_traj.xyz'.format(sys.argv[1][:-4])))
gbid        = (sys.argv[1][:-4]).split('/')[-1]
strain_mask = [0,0,1,0,0,0]
ucf         = UnitCellFilter(grain, strain_mask)
opt         = FIRE(grain)
for i in range(32):
  opt.run(fmax=0.008, steps=32)
  out.write(grain)
  if max(np.sum(grain.get_forces()**2, axis=1)**0.5) < 0.008:
    break
out.close()

E_gb = grain.get_potential_energy()
cell = grain.get_cell()
A    = cell[0][0]*cell[1][1]
H    = cell[2][2]
#Calculation dumps total energyenergy and grainboundary area data to json file.
gb_dict = {'gbid':gbid, 'E_gb':E_gb, 'E_gb_init':E_gb_init, 'A': A, 'H':H, 'n_at':len(grain), 
           'param_file':pot_file}
#add keys  
with open('subgb.json', 'w') as outfile:
  for key, value in gb_dict.items():
    j_dict[key] = value
  json.dump(j_dict, outfile, indent=2)

