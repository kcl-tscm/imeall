import argparse
import numpy as np
import operator
import sys
import transformations as quat

from fractions import gcd

class QuaternionGB(object):
  """:class:`imeall.slabmaker.gengb_from_quat.QuaternionGB` contains routines to generate 
  lists of valid grain boundary orientations (i.e. the boundary planes of
  grains for different orientation axes and misorientation angles) using
  quaternion algebra. Currently supports routines for symmetric tilt 
  grain boundaries and twist boundaries.
  """

  def __init__(self):
    pass

  def or_axis_to_quat(self, orientation_axis):
    """Given an orientation axis converts the vector to quaternion
    representation defining the orientation axis and a second quaternion
    giving the 0 of the angle. :
    
    Args:
      orientation_axis(list): orientation axis.

    Returns:
      quaternion, quaternion: planequat_1, planequat_2
    """
    #001 axis
    if (np.array(orientation_axis) == np.array([0.,0.,1])).all():
      planequat_1      = np.array([0,0,1,0])
      planequat_2      = np.array([0,1,0,0])
    #110 axis
    elif (np.array(orientation_axis) == np.array([1.,1.,0])).all():
      planequat_1      = np.array([0,0,0,1])
      planequat_2      = np.array([0,1,-1,-1])
    #111 axis
    elif (np.array(orientation_axis) == np.array([1.,1.,1])).all():
      planequat_1      = np.array([0,1,1,-2])
      planequat_2      = np.array([0,1, -1,0])
    #112 axis
    elif (np.array(orientation_axis) == np.array([1.,1.,-2])).all():
      planequat_1      = np.array([0,1,-1,0])
      planequat_2      = np.array([0,1,1,1])
    else:
      sys.exit('Orientation axis not supported please add defining quaternions.')
    return planequat_1, planequat_2

  def or_axis_to_angle(self, m, n, orientation_axis):
    """To generate an approximate spanning set of angles iterate
    over integers m and n and determine the angle from a combination
    of these integers according to the chosen orientation axis.

    Args:
      m(int): Integer rotation 
      n(int): Integer rotation

    Returns:
      float: misorientation angle from integers and orientation axis.
    """

    if (np.array(orientation_axis) == np.array([0.,0.,1])).all():
      angle = np.arccos((m**2-n**2)/(m**2+n**2))
    elif (np.array(orientation_axis) == np.array([1.,1.,0])).all():
      angle = np.arccos((m**2-2*n**2)/(m**2+2*n**2))
    elif (np.array(orientation_axis) == np.array([1.,1.,1])).all():
      angle = np.arccos((m**2-3*n**2)/(m**2+3*n**2))
    elif (np.array(orientation_axis) == np.array([1.,1.,-2])).all():
      angle = np.arccos((m**2-n**2 -(2*n)**2)/(m**2+n**2+(2*n)**2))
    else:
      sys.exit()
    return angle

  def gen_sym_tilt(self, orientation_axis=[0.,0.,1.]):
    """Generate symmetric tilt grain boundary list.
    Generate a list of the misorientation angle for the given orientation axis
    and the relevant boundary plane to generate the symmetric tilt grain
    boundary bicrystal with :func:`imeall.slabmaker.build_tilt_sym_gb`.

    Args: 
      orientation_axis(list[int]): 3d array of floats defining orientation axis.

    Returns:
      dict: Dictionary of 'x' axis co-ordinates keyed by angle rounded to two floating point numbers.
    """
    planequat_1, planequat_2 = self.or_axis_to_quat(orientation_axis)
    deg_list = []
    deg_list_py_tmp = []
    #Generate a list of (angle, boundary plane) pairs for the chosen orientation axis.
    for n  in np.arange(1,20):
      for m in np.arange(1,20):
        n = float(n)
        m = float(m)
        angle = self.or_axis_to_angle(m, n, orientation_axis)
        #The quaternion characterizing the rotation matrix:
        rotquat = quat.quaternion_about_axis(angle, orientation_axis)
        rotmquat = (1./(np.array([a for a in rotquat if a!=0.]).min())*rotquat).round(8)
        n1 = quat.quaternion_multiply(planequat_1, rotmquat)
        n2 = quat.quaternion_multiply(planequat_2, rotmquat)
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
        grcd = reduce(gcd, abs_bp1)
        bp1 = bp1/grcd
        deg = '\t [np.pi*({0}/180.), np.array([{1}, {2}, {3}])],'.format(round(180.*angle/np.pi,2), bp1[1].round(2), 
                                                                               bp1[2].round(2), bp1[3].round(2))
        if all(bp1 < 1e3):#handles floating point errors
          deg_list.append([round(180.*angle/np.pi, 2),  [deg]])
          deg_list_py_tmp.append([round(180.*angle/np.pi,2), np.array([bp1[1].round(2), bp1[2].round(2), bp1[3].round(2)])])
    deg_list = sorted(deg_list, key=lambda x: x[0])
    deg_list_py_tmp = sorted(deg_list_py_tmp, key=lambda x: x[0])
    deg_redun =  [] 
    deg_dict  = {}
    deg_list_py = []
    for deg, deg_py in zip(deg_list, deg_list_py_tmp):
      try:
        x = deg_dict[round(deg[0],2)]
      except KeyError:
       deg_dict[round(deg[0],2)] = deg[1][0]
       deg_list_py.append(deg_py)
    return deg_list_py

  def gen_sym_twist(self, orientation_axis=[0,0,1]):
    """Generate symmetric twist grain boundary list.
    Generate a list of the misorientation angle for the given orientation axis
    and the relevant boundary plane to generate the symmetric twist grain
    boundary bicrystal with :func:`imeall.slabmaker.build_twist_sym_gb`.

    Args:
      orientation_axis(list[ints]): orientation axis in ints.

    Returns:
      dict: Dictionary of 'x' axis co-ordinates keyed by angle rounded to two floating point numbers.
    """
    planequat_1, planequat_2 = self.or_axis_to_quat(orientation_axis)
    deg_list = []
    deg_list_py_tmp = []
    #Generate a list of (angle, rotated plane) pairs for the chosen orientation axis.
    for n  in np.arange(1,20):
      for m in np.arange(1,20):
        n = float(n)
        m = float(m)
        angle = self.or_axis_to_angle(m, n, orientation_axis)
        #The quaternion characterizing the rotation matrix:
        rotquat = quat.quaternion_about_axis(angle, orientation_axis)
        rotmquat = (1./(np.array([a for a in rotquat if a!=0.]).min())*rotquat).round(8)
        n1 = quat.quaternion_multiply(planequat_1, rotmquat)
        n2 = quat.quaternion_multiply(planequat_2, rotmquat)
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
        grcd = reduce(gcd, abs_bp1)
        bp1 = bp1/grcd
        deg = '\t [np.pi*({0}/180.), np.array([{1}, {2}, {3}])],'.format(round(180.*angle/np.pi,2), bp1[1].round(2), 
                                                                               bp1[2].round(2), bp1[3].round(2))
        if all(bp1 < 1e3):#handles floating point errors
          deg_list.append([round(180.*angle/np.pi, 2),  [deg]])
          deg_list_py_tmp.append([round(180.*angle/np.pi,2), np.array([bp1[1].round(2), bp1[2].round(2), bp1[3].round(2)])])
    deg_list = sorted(deg_list, key=lambda x: x[0])
    deg_list_py_tmp = sorted(deg_list_py_tmp, key=lambda x: x[0])
    deg_redun =  [] 
    deg_dict  = {}
    deg_list_py = []
  #only print unique boundaries
    for deg, deg_py in zip(deg_list, deg_list_py_tmp):
      try:
        x = deg_dict[round(deg[0], 2)]
      except KeyError:
       deg_dict[round(deg[0],2)]=deg[1][0]
       deg_list_py.append(deg_py)
    return deg_list_py

if __name__=='__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--gb_type','-gbt', help='grain boundary type (tilt, twist).', default='tilt')
  parser.add_argument('--or_axis','-or', help='orientation axis', nargs='+', default=[1,1,0])
  args = parser.parse_args()
  or_axis = map(int, args.or_axis)
  quatgb = QuaternionGB()

  if args.gb_type=='tilt':
    quatgb.gen_sym_tilt(orientation_axis=or_axis)
  elif args.gb_type=='twist':
    quatgb.gen_sym_twist(orientation_axis=or_axis)
  else:
    print 'grain boundary type not recognized.'

