import os
import sys
import json
import shutil
import ase.io        
import numpy as np
from pprint import pprint
from cStringIO           import StringIO
from ase.optimize.sciopt import SciPyFminCG
from quippy              import Atoms, Potential, frange
from ase.constraints     import UnitCellFilter, StrainFilter
from quippy.io           import AtomsWriter, AtomsReader, write
from ase.optimize        import BFGS, FIRE, LBFGS, MDMin, QuasiNewton

def relax_gb(gb_file='file.xyz'):
  def converged(grain, smax, fmax):
    maxstress = max(grain.get_stress().ravel())
    rmsforces = np.sum(grain.get_forces()**2, axis=1)**0.5
    maxforce  = max(rmsforces)
    if maxforce < fmax and maxstress < smax:
      return True
    return False
  
  with open('subgb.json', 'r') as outfile:
    j_dict = json.load(outfile)
  
  #Rescaled to 2.83 which is magic number for the grain boundary canonical case.
  #Ada
  try:
    #POT_DIR = '/users/k1511981/pymodules/imeall/imeall/potentials' 
    POT_DIR     = os.environ['POTDIR']
  except:
    sys.exit("PLEASE SET export POTDIR='path/to/potfiles/'")
  #Retina
  #POT_DIR = '/Users/lambert/pymodules/imeall/imeall/potentials' 
  try: 
    param_file = j_dict['param_file']
    if param_file == 'iron_mish.xml':
      eam_pot = os.path.join(POT_DIR, 'iron_mish.xml')
      r_scale = 1.0129007626
    elif param_file == 'Fe_Mendelev.xml':
      eam_pot = os.path.join(POT_DIR, 'Fe_Mendelev.xml')
      r_scale = 1.00894848312
    elif param_file == 'PotBH.xml':
      eam_pot = os.path.join(POT_DIR, 'PotBH.xml')
      r_scale = 1.00894848312
    elif param_file == 'Fe_Ackland.xml':
      eam_pot = os.path.join(POT_DIR,'Fe_Ackland.xml')
      r_scale = 1.00894185389
    elif param_file == 'Fe_Dudarev.xml':
      eam_pot = os.path.join(POT_DIR,'Fe_Dudarev.xml')
      r_scale = 1.01279093417 
  except KeyError:
    print 'No EAM pot relax failed!'
    sys.exit()
  
  print 'Using: ', eam_pot
  pot_file    = eam_pot.split('/')[-1]
  #grain       = Atoms('{0}'.format(sys.argv[1]))
  print '{0}.xyz'.format(gb_file)
  print os.getcwd()
  grain       = Atoms('{0}.xyz'.format(gb_file))
  pot         = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)
  grain.set_calculator(pot)
  E_gb_init   = grain.get_potential_energy()
  alpha       = E_gb_init
  traj_file   = gb_file
  out         = AtomsWriter('{0}'.format('{0}_traj.xyz'.format(traj_file)))
  #gbid        = (gb_file[:-4]).split('/')[-1]
  gbid        = gb_file
  strain_mask = [0,0,1,0,0,0]
  ucf         = UnitCellFilter(grain, strain_mask)
  opt         = FIRE(ucf)

  cell = grain.get_cell()
  A    = cell[0][0]*cell[1][1]
  H    = cell[2][2]
  #Calculation dumps total energyenergy and grainboundary area data to json file.
  gb_dict = {'gbid':gbid, 'E_gb_init':E_gb_init, 'A': A, 'H':H, 'n_at':len(grain), 
             'param_file':pot_file, 'converged': False}
  #Write an initial dict so we know if the system has been initialized but the calculation is not finished.
  with open('subgb.json', 'w') as outfile:
    for key, value in gb_dict.items():
      j_dict[key] = value
    json.dump(j_dict, outfile, indent=2)

  for i in range(32):
    opt.run(fmax=0.01, steps=200)
    out.write(grain)
    if max(np.sum(grain.get_forces()**2, axis=1)**0.5) < 0.01:
      CONVERGED = True
      break

  out.close()
  E_gb = grain.get_potential_energy()
  gb_dict = {'gbid':gbid, 'E_gb':E_gb, 'E_gb_init':E_gb_init, 'A': A, 'H':H, 'n_at':len(grain), 
             'param_file':pot_file, 'converged':CONVERGED}
  #Add keys  
  with open('subgb.json', 'w') as outfile:
    for key, value in gb_dict.items():
      j_dict[key] = value
    json.dump(j_dict, outfile, indent=2)

if __name__ == '__main__':
  relax_gb()
