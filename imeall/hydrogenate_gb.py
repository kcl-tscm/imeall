import os
import glob
import shutil
import argparse
import numpy as np

from   ase.optimize import FIRE
from   ase.optimize.precon import LBFGS
from   ase.lattice.cubic  import BodyCenteredCubic
from   ase.constraints    import UnitCellFilter, StrainFilter

from   quippy import Atoms, Potential, AtomsReader
from   fracture.hydrify_cracktips import Hydrify

j_file = '00195301120_v6bxv2_tv0.4bxv0.3_d1.9z_traj.xyz'
E_H2   = -4.73831215118
E_diss =  0.291385500301

hydrify = Hydrify()
gb = AtomsReader(j_file)[-1]
#gb = BodyCenteredCubic(directions = [[1,0,0], [0,1,0], [0,0,1]],
#                            size = (8,8,16), symbol='Fe', pbc=(1,1,1),
#                            latticeconstant = 2.83)
try:
  #POT_DIR = '/users/k1511981/pymodules/imeall/imeall/potentials' 
  POT_DIR     = os.environ['POTDIR']
except:
  sys.exit("PLEASE SET export POTDIR='path/to/potfiles/'")

eam_pot = os.path.join(POT_DIR, 'PotBH.xml')
r_scale = 1.00894848312
pot         = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)
cell = gb.get_cell()
A    = cell[0][0]*cell[1][1]
Height    = cell[2][2]
print 'A', A, 'H', Height

parser = argparse.ArgumentParser()
parser.add_argument("-s","--strain",  type=float, help="strain parameter", default=1.0)
parser.add_argument("-d","--dirname", help="Directory for calculation")
parser.add_argument("-b","--bulk",    help="Deposit H in bulk positions or at grainboundary.", action="store_true")
args    = parser.parse_args()
bulk    = args.bulk
strain  = args.strain
dirname = args.dirname

scratch = os.getcwd()
try:
  os.mkdir(dirname)
except:
  pass
os.chdir(dirname)

if bulk:
  z_plane = Height/2.0
else:
  z_plane = 102.119

cell       = gb.get_cell()
cell [2,2] = cell[2,2]*strain
gb.set_cell(cell, scale_atoms=True)
with open('H_en.dat', 'a') as f:
  print >> f,'n_H  d_H  E_gbh  E_gb  n_H*0.5*E(H2) E_sg   n_H*E_H  Corr_E'
  for d_H in np.arange(3.5, 10.0, 1.0):
    print 'Inter_H distance', d_H
    gb_h = gb.copy()
    if bulk:
      h_list = hydrify.hydrogenate_gb(gb, d_H=d_H, z_plane=z_plane, tetrahedral=True)
    else:
      h_list = hydrify.hydrogenate_gb(gb, d_H=d_H, z_plane=z_plane, tetrahedral=False)

    for h in h_list:
      gb_h.add_atoms(h,1)

    gb.set_calculator(pot)
    gb_h.set_calculator(pot)
    
    E_gb   = gb.get_potential_energy()
    strain_mask = [0,0,0,0,0,0]
    ucf         = UnitCellFilter(gb_h, strain_mask)
    opt         = FIRE(ucf)
    opt.run(fmax=0.05, steps=1000)
    gb_h.write('Hydrogenated_relaxed_d{0}.xyz'.format(d_H))
    E_gbh = gb_h.get_potential_energy()
    n_H  = sum([at.number == 1 for at in gb_h])
    E_sg = E_gbh-E_gb-n_H*0.5*E_H2 
    print >> f, n_H, d_H, E_gbh, E_gb, n_H*0.5*E_H2, E_sg, n_H*E_diss, E_sg-n_H*E_diss

