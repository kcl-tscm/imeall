import os
import sys
import json
from   ase.io       import Trajectory
from   ase.optimize import BFGS, FIRE
from   quippy       import Atoms

class ImeallIO(object):
# Primary IO Object for routines for searching Imeall Directory tree
# creating new grains, and subgrain directories. 
# Sub-Grain directory (each 'massaged' grain boundary will have its own
# directory) supercells of the grain count as subgrains finally a subdirectory
# of defects is included for vacancies and interstitials.
  def __init__(self):
    pass

  def make_dir(self, target_dir, dir_name):
    target_subdir = os.path.join(target_subgrain_dir)
    if not os.path.isdir(target_subgrain_dir):
      os.mkdir(target_subgrain_dir)
    else:
      print '\t directory already exists'
    return target_subdir

class GBRelax(object):
  def __init__(self, target_dir='./', gbid='0000000000', calc_type='EAM',
               param_file = 'IP EAM_ErcolAd', potential = 'iron_mish.xml'):
# Here we initialize some naming conventions, the calculation type, and
# necessary input files.
    self.gbid        =  gbid
    self.target_dir  = os.path.join(target_dir, calc_type)
    self.name        = 'subgrain'
    self.calc_type   = calc_type
    self.struct_file = ''
    self.param_file  = param_file
    self.potential   = potential
    self.fmax        = 1E-4

  def delete_atoms(self, rcut=2.0):
# Delete atoms below a certain threshold
    x = Atoms('{0}.xyz'.format(sys.argv[1]))
    x.set_cutoff(3.0)
    x.calc_connect()
    x.calc_dists()
    rem=[]
    for j in frange(x.n):
      for i in frange(j, x.n):
        if 0. < x.distance_min_image(i, j) < rcut and j!=i:
          rem.append(sorted([j,i]))
    rem = list(set([a[0] for a in rem]))
    if len(rem) >0:
      x.remove_atoms(rem)
    else:
      print 'No duplicate atoms in list.'
    self.name    = 'gbid_d{0}'.format(str(rcut)) 
    self.subgrain_dir = ImeallIO.make_dir(self.target_dir, self.name)

    self.struct_file  = os.path.join(len(rem), rcut)
    self.struct_file    = os.path.join(target_subgrain_dir,struct_file)
    x.write('{0}.xyz'.format(struct_file)

  def translate(self):
# Displace grain_a relative to grain_b
# can be accomplished by moving along x or y
# and wrapping atoms back to unit cell
    pass

# Take the grain_boundary structure and relax it
  def go_relax(self):
    if self.struct_file = '':
      print 'No pointer to struct file. Have you specified a job?'
      return None
    else:
      print 'Print Relaxing {0} with {1} potential'.format(self.struct_file, self.calc_type)
    if self.calc_type == 'EAM':
      grain = Atoms('{0}.xyz'.format())
# Load potential and Attach calculator
      pot = Potential(self.potential, parame_filename=self.param_file)
      grain.set_calculator(pot)
      strain_mask = [0,0,0,0,0,0]
      ucf         = UnitCellFilter(grain, mask=strain_mask)
      traj        = Trajectory(''.format(self.target_dir, self.name), mode='w', atoms=grain)
      opt         = BFGS(ucf)
      opt.attach(traj)
      opt.run(fmax=self.fmax)
      E_gb = grain.get_potential_energy()
      cell = grain.get_cell()
      A = cell[0][0]*cell[1][1]
      gb_dict = {'E_gb':E_gb, 'E_gb_init':E_gb_init, 'A': A}
# Calculation dumps total energyenergy and grainboundary area data.
      with(os.path.join(self.subgrain_dir, 'subgb.json')) as outfile:
        json.dumps(gb_dict, outfile)
    else:
# I suppose the idea here would be something like
# elif self.calc_type == 'DFT':
#   ssh submit xyz mira or archer or whereever 
#   sh mdir /mira/imeall/.../subgrain_dir
#   sh >>mira/imeall/.../subgrain_dir/POSCAR 
#   sh submit vasp.run
#   wait ....
#   sh get output
      pass

