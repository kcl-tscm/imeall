from __future__ import print_function
import os
import sys

import json 
import numpy as np

from ase.optimize import FIRE
from ase.optimize.precon import PreconFIRE, Exp
import argparse

from distutils import spawn

from imeall import app
from ForceMixerCarver import ForceMixingCarvingCalculator
from matscipy.socketcalc import VaspClient, SocketCalculator

from quippy import AtomsReader, AtomsWriter
from quippy import Potential

mpirun = spawn.find_executable('mpirun')
#vasp = '/home/eng/essswb/vasp.5.4.1.patched/bin/vasp_std'
vasp = '/home/mmm0007/vasp/vasp.5.4.1/bin/vasp_std'

parser = argparse.ArgumentParser()
parser.add_argument("--buff", "-b", type=float, default=6.0)
parser.add_argument("--qm_radius", "-q", type=float, default=2.5)
parser.add_argument("--use_socket", "-u", type=bool, default=True)
args=parser.parse_args()

buff = args.buff
qm_radius = args.qm_radius

POT_DIR = os.path.join(app.root_path, 'potentials')
eam_pot = os.path.join(POT_DIR, 'PotBH.xml')
r_scale = 1.00894848312
mm_pot = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)
#r_scale = 1.0089
#mm_pot_mod = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)



with open('subgb.json','r') as f:
  subgb_dict = json.load(f)
#[11.79728921, 4.61819134, 34.51540462]
h_pos = subgb_dict['site']
print (h_pos)
#disloc_fin, disloc_ini, bulk = sd.make_barrier_configurations(calculator=lammps, cylinder_r=cylinder_r)
gb_cell = AtomsReader('interstitial_traj.xyz')[-1]

vasp_args = dict(xc='PBE', amix=0.01, amin=0.001, bmix=0.001, amix_mag=0.01, bmix_mag=0.001,
                 kpts=[1, 1, 1], kpar=1, lreal='auto', nelmdl=-15, ispin=2, prec='Accurate',
                 nelm=100, algo='VeryFast', lplane=False, lwave=False, lcharg=True, istart=0, 
                 magmom=2.6*len(gb_cell), maxmix=30, #https://www.vasp.at/vasp-workshop/slides/handsonIV.pdf #for badly behaved clusters.
                 voskown=0, ismear=1, sigma=0.1, isym=2) # possibly try iwavpr=12, should be faster if it works

x, y, z = gb_cell.positions.T
radius1 = np.sqrt((x - h_pos[0])**2 + (y-h_pos[1])**2 + (z-h_pos[2])**2)

qm_region_mask = (radius1 < qm_radius)
qm_buffer_mask = (radius1 < qm_radius + buff)

print ("\nNumber of atoms in qm region of %.1f" % qm_radius +
                                "A : %i" % np.count_nonzero(qm_region_mask))

print ("together with the buffer of %.1f" % (qm_radius + buff ) +
                                "A %i" % np.count_nonzero(qm_buffer_mask))


#parallel config.
procs = 24
kpts = [1, 1, 1]
# need to have procs % n_par == 0
n_par = 1
if procs <= 8:
  n_par = procs
else:
  for _npar in range(2, int(np.sqrt(1.*procs))):
    if procs % int(_npar) == 0:
      n_par = procs // int(_npar)


if args.use_socket:
  vasp_client = VaspClient(client_id=0, npj=procs, ppn=1,
                         exe=vasp, mpirun=mpirun, parmode='mpi',
                         ibrion=13, nsw=1000000,
                         npar=n_par, **vasp_args)

  qm_pot = SocketCalculator(vasp_client)

else:
  pass

qmmm_pot = ForceMixingCarvingCalculator(gb_cell, qm_region_mask,
                                        mm_pot, #mm_pot_mod, #for testing
                                        qm_pot, 
                                        buffer_width=buff,
                                        pbc_type=[False, False, False])

#opt = PreconFIRE(gb_cell, precon=Exp())
def pass_trajectory_context(trajectory, dynamics):
  def traj_writer(dynamics):
      trajectory.write(dynamics.atoms)
  return traj_writer
trajectory = AtomsWriter('gb_traj.xyz')

opt = FIRE(gb_cell)
opt.attach(pass_trajectory_context(trajectory, opt),1,opt)
gb_cell.set_calculator(qmmm_pot)
opt.run(fmax=1.0e-3)
gb_cell.write('dft_relaxed.xyz')
sock_calc.shutdown()

