import os
import sys
import json
import shutil
from   cStringIO       import StringIO
import ase.io        
from   ase.constraints import UnitCellFilter, StrainFilter
from   ase.optimize    import BFGS, FIRE, LBFGS, MDMin, QuasiNewton
from   ase.optimize.sciopt import SciPyFminCG
from   quippy              import Atoms, Potential, frange
from   quippy.io           import AtomsWriter, AtomsReader, write
import numpy as np

from pprint import pprint


def converged(grain, smax, fmax):
  maxstress = max(grain.get_stress().ravel())
  rmsforces = np.sum(grain.get_forces()**2, axis=1)**0.5
  maxforce  = max(rmsforces)
  if maxforce < fmax and maxstress < smax:
    return True
  return False

eam_pot = '/users/k1511981/sharedscratch/grain_boundaries/iron_mish.xml'
grain       = Atoms('{0}'.format(sys.argv[1]))
pot         = Potential('IP EAM_ErcolAd', param_filename=eam_pot)
grain.set_calculator(pot)
E_gb_init   = grain.get_potential_energy()

alpha       = E_gb_init
out         = AtomsWriter('{0}'.format('{0}_traj.xyz'.format(sys.argv[1][:-4])))
gbid        = (sys.argv[1][:-4]).split('/')[-1]
#mdmin_cell  = MDMin(strain_mask,dt=0.1)
strain_mask = [0,0,1,0,0,0]
ucf         = UnitCellFilter(grain)
opt         = FIRE(grain)
#print grain.get_stress()
#while not converged(grain, smax=0.05, fmax=0.05):
for i in range(32):
  opt.run(fmax=0.01, steps=25)
  out.write(grain)
  if max(np.sum(grain.get_forces()**2,axis=1)**0.5) < 0.01:
    break

out.close()
E_gb = grain.get_potential_energy()
cell = grain.get_cell()
A    = cell[0][0]*cell[1][1]
H    = cell[2][2]
# Calculation dumps total energyenergy and grainboundary area data to json file.
gb_dict = {'gbid':gbid,'E_gb':E_gb, 'E_gb_init':E_gb_init, 'A': A, 'H':H, 'n_at':len(grain)}
#add keys  
with open('subgb.json', 'r') as outfile:
  j_dict = json.load(outfile)
with open('subgb.json', 'w') as outfile:
  for key, value in gb_dict.items():
    j_dict[key] = value
  json.dump(j_dict, outfile, indent=2)
