from ase.io import read
from ase.io import write
from ase.lattice.spacegroup import crystal
from ase.lattice.surface import surface, bcc111,bcc110
from ase.lattice.cubic   import BodyCenteredCubic
from ase.utils.geometry  import get_duplicate_atoms

import os
import json
import time
import operator
import numpy as np
import transformations as quat

from slabmaker import gen_csl
from fractions import gcd
from quippy import io
from quippy import set_fortran_indexing
from quippy import Atoms
try:
  from qlab import set_fortran_indexing, view, gcat
except:
  pass

set_fortran_indexing(False)

def integralize_quaternion(n1, n2, angle):
  """
  :method: Integralize quaternion.
  :return: formatted string of angle and appropriate x vector.
  """
  comm_denom = []
  for a in n1:
    if (a%1).round(1) != 0:
      comm_denom.append(1./(np.abs(a)%1))
  comm_denom = np.array(comm_denom)
  if len(comm_denom)!=0:
    bp1 = (comm_denom.max()*n1).round(2)
    bp2 = (comm_denom.max()*n2).round(2)
  else:
    bp1 = n1
    bp2 = n2
  abs_bp1 = map(np.abs, bp1)
  grcd  = reduce(gcd, abs_bp1)
  bp1     = bp1/grcd
  rad  = '\t [{0}, np.array([{1} {2} {3}])],'.format(round(angle, 2), bp1[1].round(0), bp1[2].round(0), bp1[3].round(0))
  deg  = '\t [np.pi*({0}/180.), np.array([{1}, {2}, {3}])],'.format(round(180.*angle/np.pi,2), bp1[1].round(2), bp1[2].round(2), bp1[3].round(2))
  return [round(180.*angle/np.pi, 2),  [deg]]

def gen_sym_twist_gb(or_axis=[0,0,1]):
  """
  :method: initialize symmetric twist grain boundary
  """
  if (np.array(or_axis) == np.array([0,0,1])).all():
    planequat_1      = np.array([0,0,1,0])
    planequat_2      = np.array([0,1,0,0])
    #planequat_1      = np.array([0,1,0,0])
    #planequat_2      = np.array([0,0,1,0])
  print planequat_1
  print planequat_2
  deg_list = []
  for n  in np.arange(1,20):
    for m in np.arange(1,20):
      print n,m
      n = float(n)
      m = float(m)
      if (np.array(or_axis) == np.array([0,0,1])).all() :
        angle = np.arccos((m**2-n**2)/(m**2+n**2))
      else:
        print '{} not implemented'.format(or_axis)
        sys.exit()
      rotm     = quat.rotation_matrix(angle, map(float,or_axis))
      rotquat  = quat.quaternion_from_matrix(rotm).round(8)
      rotmquat = (1./(np.array([a for a in rotquat if a!=0.]).min())*rotquat).round(5)
      n1 = quat.quaternion_multiply(planequat_1, rotmquat)
      n2 = quat.quaternion_multiply(planequat_2, rotmquat)
      deg_list.append(integralize_quaternion(n1, n2, angle))
  deg_list = sorted(deg_list, key=lambda x: x[0])
  deg_dict = {}
#only print unique boundaries
  for deg in deg_list:
    try:
      x = deg_dict[round(deg[0],2)]
    except KeyError:
     deg_dict[round(deg[0],2)]=deg[1][0]
     print deg[1][0]

def build_twist_sym_gb(gbid='', bp =[0,0,1], v=[3,5,0],
                      c_space=None, target_dir=None, rbt = None):
  """
  :method:`build_twist_sym_gb` the or_axis is the boundary plane in this case.
  the X axis of the grain.
  """
  print v
  bpxv    = [(bp[1]*v[2]-v[1]*bp[2]), (bp[2]*v[0]-bp[0]*v[2]), (bp[0]*v[1]- v[0]*bp[1])]
  grain_a =  BodyCenteredCubic(directions = [v, bpxv, bp],
                               size = (1,1,1), symbol='Fe', pbc=(1,1,1),
                               latticeconstant = 2.83)
  n_grain_unit = len(grain_a)
  n = 2
  while(grain_a.get_cell()[2,2]< 11.0 and n < 10):
    grain_a = BodyCenteredCubic(directions = [v, bpxv, bp],
                              size = (1,1,n), symbol='Fe', pbc=(1,1,1),
                              latticeconstant = 2.83)
    n += 1
  v2      = v.copy()
  v2[0]   = -v2[0]
  print v2
  bpxv    = [(bp[1]*v2[2]-v2[1]*bp[2]),(bp[2]*v2[0]-bp[0]*v2[2]),(bp[0]*v2[1]- v2[0]*bp[1])]
  grain_b = BodyCenteredCubic(directions = [v2, bpxv, bp],
                              size = (1,1,1), symbol='Fe', pbc=(1,1,1),
                              latticeconstant = 2.83)
  n = 2
  while(grain_b.get_cell()[2,2]< 22.0 and n < 10):
    grain_b = BodyCenteredCubic(directions = [v2, bpxv, bp],
                              size = (1,1,n), symbol='Fe', pbc=(1,1,1),
                              latticeconstant = 2.83)
    n += 1
  grain_c = grain_a.copy()
  if c_space==None:
    s1      = surface('Fe', (0,0,1), n)
    c_space = grain_b.get_cell()[2,2]/float(n)
    c_space = np.array(2.83)
    s2 = surface('Fe', (map(int, v)), 1)
    x_space = s2.get_cell()[0,0] #-s1.positions[:,2].max()
    s3 = surface('Fe', (map(int, bpxv)), 1)
    y_space = s3.get_cell()[1,1] #-s1.positions[:,2].max()
  print '\t Interplanar spacing: ', x_space.round(2), y_space.round(2), c_space.round(2), 'A'
  if sum([int(n)%2 for n in bp])%2 == 0 :
    grain_b.positions[:,2]  -= (grain_b.positions[:,2].max() + c_space)
  else:
    grain_b.positions[:,2]  -= (grain_b.positions[:,2].max() + c_space/2.0)
  grain_c.extend(grain_b)
  grain_c.set_cell([grain_c.get_cell()[0,0], grain_c.get_cell()[1,1], 2.*grain_c.get_cell()[2,2]])
  grain_c.positions[:,2] += abs(grain_c.positions[:,2].min())
# In BCC lattice spacing is different depending on whether the
# z-plane has an even or odd number of odd miller indices:
  if sum([int(n)%2 for n in bp])%2 == 0 :
    grain_a.positions[:,2] -= (grain_a.positions[:,2].max()+c_space)
  else:
    grain_a.positions[:,2] -= (grain_a.positions[:,2].max()+c_space/2.0)
  grain_c.extend(grain_a)
#grain_c.set_cell([grain_c.get_cell()[0,0], grain_c.get_cell()[1,1], 2.*grain_c.get_cell()[2,2]])
  grain_c.positions[:,2] += abs(grain_c.positions[:,2].min())
  dups = get_duplicate_atoms(grain_c)
  for dup in dups:
      grain_c[dup[1]].position[2] += grain_a.get_cell()[2,2]
  shared_atoms = [grain_c[dup[0]] for dup in dups]
  for at in shared_atoms:
    print '\t', at.index, at.position
  print '\t There are {0} shared atoms'.format(len(shared_atoms))
#Space bicrystal properly.
  if sum([int(n)%2 for n in bp])%2 == 0 :
	  grain_c.center(vacuum=c_space/2.0, axis=2)
  else:
    grain_c.center(vacuum=c_space/4.0, axis=2)
  grain_c.info['adsorbate_info']=None
  grain_c = Atoms(grain_c)
  if target_dir != None:
    print '\t Writing {0}.xyz to file'.format(gbid)
    n_grain_unt = len(grain_c)
    z_planes = [round(atom.position[2],4) for atom in shared_atoms]
    z_planes = list(sorted(set(z_planes)))
    return [z_planes, len(dups), n_grain_unit, grain_c]
  else:
    return grain_c

if __name__=='__main__':
  sym_twist_001 = [
	         [np.pi*(7.15/180.), np.array([1.0, 16.0, 0.0])],
	         [np.pi*(7.63/180.), np.array([1.0, 15.0, 0.0])],
					 [np.pi*(8.17/180.), np.array([1.0, 14.0, 0.0])],
					 [np.pi*(8.8/180.), np.array([1.0, 13.0, 0.0])],
					 [np.pi*(9.53/180.), np.array([1.0, 12.0, 0.0])],
					 [np.pi*(10.39/180.), np.array([1.0, 11.0, 0.0])],
					 [np.pi*(11.42/180.), np.array([1.0, 10.0, 0.0])],
					 [np.pi*(12.68/180.), np.array([1.0, 9.0, 0.0])],
					 [np.pi*(14.25/180.), np.array([1.0, 8.0, 0.0])],
					 [np.pi*(15.19/180.), np.array([2.0, 15.0, 0.0])],
					 [np.pi*(16.26/180.), np.array([1.0, 7.0, 0.0])],
					 [np.pi*(17.49/180.), np.array([2.0, 13.0, 0.0])],
					 [np.pi*(18.92/180.), np.array([1.0, 6.0, 0.0])],
					 [np.pi*(20.61/180.), np.array([2.0, 11.0, 0.0])],
					 [np.pi*(21.24/180.), np.array([3.0, 16.0, 0.0])],
					 [np.pi*(22.62/180.), np.array([1.0, 5.0, 0.0])],
					 [np.pi*(24.19/180.), np.array([3.0, 14.0, 0.0])],
					 [np.pi*(25.06/180.), np.array([2.0, 9.0, 0.0])],
					 [np.pi*(25.99/180.), np.array([3.0, 13.0, 0.0])],
					 [np.pi*(28.07/180.), np.array([1.0, 4.0, 0.0])],
					 [np.pi*(30.51/180.), np.array([3.0, 11.0, 0.0])],
					 [np.pi*(31.89/180.), np.array([2.0, 7.0, 0.0])],
					 [np.pi*(33.4/180.), np.array([3.0, 10.0, 0.0])],
					 [np.pi*(34.21/180.), np.array([4.0, 13.0, 0.0])],
					 [np.pi*(34.71/180.), np.array([5.0, 16.0, 0.0])],
					 [np.pi*(36.87/180.), np.array([1.0, 3.0, 0.0])],
					 [np.pi*(39.31/180.), np.array([5.0, 14.0, 0.0])],
					 [np.pi*(41.11/180.), np.array([3.0, 8.0, 0.0])],
					 [np.pi*(43.6/180.), np.array([2.0, 5.0, 0.0])],
					 [np.pi*(45.24/180.), np.array([5.0, 12.0, 0.0])],
					 [np.pi*(46.4/180.), np.array([3.0, 7.0, 0.0])],
					 [np.pi*(47.26/180.), np.array([7.0, 16.0, 0.0])],
					 [np.pi*(47.92/180.), np.array([4.0, 9.0, 0.0])],
					 [np.pi*(48.89/180.), np.array([5.0, 11.0, 0.0])],
					 [np.pi*(49.55/180.), np.array([6.0, 13.0, 0.0])],
					 [np.pi*(50.03/180.), np.array([7.0, 15.0, 0.0])],
					 [np.pi*(53.13/180.), np.array([1.0, 2.0, 0.0])],
					 [np.pi*(58.11/180.), np.array([5.0, 9.0, 0.0])],
					 [np.pi*(61.93/180.), np.array([3.0, 5.0, 0.0])],
					 [np.pi*(64.94/180.), np.array([7.0, 11.0, 0.0])],
					 [np.pi*(67.38/180.), np.array([2.0, 3.0, 0.0])],
					 [np.pi*(69.39/180.), np.array([9.0, 13.0, 0.0])],
					 [np.pi*(71.08/180.), np.array([5.0, 7.0, 0.0])],
					 [np.pi*(72.51/180.), np.array([11.0, 15.0, 0.0])],
					 [np.pi*(73.74/180.), np.array([3.0, 4.0, 0.0])],
					 [np.pi*(75.75/180.), np.array([7.0, 9.0, 0.0])],
					 [np.pi*(77.32/180.), np.array([4.0, 5.0, 0.0])],
					 [np.pi*(78.58/180.), np.array([9.0, 11.0, 0.0])],
					 [np.pi*(79.61/180.), np.array([5.0, 6.0, 0.0])],
					 [np.pi*(80.47/180.), np.array([11.0, 13.0, 0.0])],
					 [np.pi*(81.2/180.), np.array([6.0, 7.0, 0.0])],
					 [np.pi*(81.83/180.), np.array([13.0, 15.0, 0.0])],
					 [np.pi*(82.37/180.), np.array([7.0, 8.0, 0.0])],
					 [np.pi*(83.27/180.), np.array([8.0, 9.0, 0.0])],
					 [np.pi*(83.97/180.), np.array([9.0, 10.0, 0.0])],
					 [np.pi*(84.55/180.), np.array([10.0, 11.0, 0.0])],
					 [np.pi*(85.02/180.), np.array([11.0, 12.0, 0.0])],
					 [np.pi*(85.42/180.), np.array([12.0, 13.0, 0.0])],
					 [np.pi*(85.76/180.), np.array([13.0, 14.0, 0.0])],
					 [np.pi*(86.05/180.), np.array([14.0, 15.0, 0.0])],
					 [np.pi*(86.3/180.), np.array([15.0, 16.0, 0.0])],
	         [np.pi*(90.0/180.), np.array([1.0, 1.0, 0.0])]]

  #gen_sym_twist_gb()
  #build_twist_sym_gb(target_dir='./')
  orientation_axis  = np.array([0,0,1])
  #surfaces = [[np.pi*(0.0), np.array([0,0,1])]]
  #surfaces = [[np.pi*(61.93/180.), np.array([3.0, 5.0, 0.0])]]
  surfaces = sym_twist_001
  for gb in surfaces:
    angle_str      = str(round((gb[0]*180./np.pi),2)).replace('.', '')
    if len(angle_str) > 4:
      angle_str = angle_str[:-1]
    elif len(angle_str) < 4:
      angle_str = angle_str + '0'
    gbid = '{0}{1}{2}'.format(orientation_axis[0], orientation_axis[1],orientation_axis[2]) \
           + angle_str  \
           + '{0}{1}{2}'.format(orientation_axis[0],orientation_axis[1], orientation_axis[2])
    print '\t Grain Boundary ID',  gbid
    gb_dir     = os.path.join('./boundaries/twist', '001')
    target_dir = os.path.join(gb_dir, gbid)
    print '\t Grain Boundary Dir', gb_dir
    if not os.path.isdir(target_dir):
      os.mkdir(target_dir)
    else:
      'directory already exists'
    gen_csl(orientation_axis, gb, target_dir=target_dir, gbid=gbid, gb_type="twist")
    zplanes, dups, nunitcell, grain_c = build_twist_sym_gb(gbid, bp=orientation_axis, v=gb[1], 
                                                           c_space=None, target_dir=target_dir,
                                                           rbt=[0.0, 0.0])
    cell    = grain_c.get_cell()
    A       = cell[0][0]*cell[1][1]
    H       = cell[2][2]
    gb_dict = {"gbid":gbid, "boundary_plane":list(orientation_axis),
               "orientation_axis":list(orientation_axis), 
               "type": "symmetric twist boundary",
               "angle": gb[0], "zplanes":zplanes, "coincident_sites": dups,
               "n_at": nunitcell, 'A':A , 'H':H}

    with open(os.path.join(target_dir, 'gb.json'), 'w') as outfile:
      json.dump(gb_dict, outfile, indent=2)

    grain_c.write(os.path.join(target_dir, '{}.xyz'.format(gbid)))

