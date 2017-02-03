import sys
import numpy as np
import transformations as quat
import operator
from   fractions import gcd

class QuaternionGB(object):
  """
  This object contains routines to generate lists of valid 
  grain boundary orientations (i.e. the boundary planes of 
  grain for different orientation axes and misorientation angles) using
  quaternion algebra. Routines for symmetric tilt, and twist grain
  boundaries.
  """
  def __init__(self):
    pass

  def gen_sym_tilt(self, orientation_axis=[0.,0.,1.]):
    """
    Generate symmetric tilt grain boundary list.
    :attributes: 
      orientation_axis: 3d array of floats.
    :returns: A printed list of the misorientation 
              angle for the given orientation axis
              and the relevant tilt boundary plane so the grain
              boundary bicrystal can be generated with 
              slabmaker.build_tilt_sym_gb.
    """
    #001 axis
    if (np.array(orientation_axis) == np.array([0.,0.,1])).all():
      planequat_1      = np.array([0,0,1,0])
      planequat_2      = np.array([0,1,0,0])
    #110 axis
    elif (np.array(orientation_axis) == np.array([1.,1.,0])).all():
      planequat_1      = np.array([0,0,0,1])
      planequat_2      = np.array([0,2,-2,-2])
    #111 axis
    elif (np.array(orientation_axis) == np.array([1.,1.,1])).all():
      planequat_1      = np.array([0,1,1,-2])
      planequat_2      = np.array([0,1, -1,0])
    elif (np.array(orientation_axis) == np.array([1.,1.,-2])).all():
      planequat_1      = np.array([0,1,-1,0])
      planequat_2      = np.array([0,1,1,1])
    
    print planequat_1
    print planequat_2
    
    deg_list = []
#Generate a list of (angle, boundary plane) pairs for
#the chosen orientation axis.
    for n  in np.arange(1,20):
      for m in np.arange(1,20):
        n = float(n)
        m = float(m)
        if (np.array(orientation_axis) == np.array([0.,0.,1])).all() :
          angle = np.arccos((m**2-n**2)/(m**2+n**2))
        elif (np.array(orientation_axis) == np.array([1.,1.,0])).all() :
          angle = np.arccos((m**2-2*n**2)/(m**2+2*n**2))
        elif (np.array(orientation_axis) == np.array([1.,1.,1])).all() :
          angle = np.arccos((m**2-3*n**2)/(m**2+3*n**2))
        elif (np.array(orientation_axis) == np.array([1.,1.,-2])).all() :
          angle = np.arccos((m**2-n**2 -(2*n)**2)/(m**2+n**2+(2*n)**2))
    #The quaternion characterizing the rotation matrix:
        rotm     = quat.rotation_matrix(angle, orientation_axis)
        rotquat  = quat.quaternion_from_matrix(rotm).round(8)
        rotmquat = (1./(np.array([a for a in rotquat if a!=0.]).min())*rotquat).round(5)
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
        grcd  = reduce(gcd, abs_bp1)
        bp1     = bp1/grcd
        rad  = '\t [{0}, np.array([{1} {2} {3}])],'.format(round(angle, 2),bp1[1].round(0), bp1[2].round(0), bp1[3].round(0))
        deg  = '\t [np.pi*({0}/180.), np.array([{1}, {2}, {3}])],'.format(round(180.*angle/np.pi,2), bp1[1].round(2), bp1[2].round(2), bp1[3].round(2))
        deg_list.append([round(180.*angle/np.pi, 2),  [deg]])
    deg_list = sorted(deg_list, key=lambda x: x[0])
    deg_redun =  [] 
    deg_dict  = {}
#print unique angles:
    for deg in deg_list:
      try:
        x=deg_dict[round(deg[0],2)]
      except KeyError:
       deg_dict[round(deg[0],2)]=deg[1][0]
       print deg[1][0]

if __name__=='__main__':
 quatgb = QuaternionGB() 
 quatgb.gen_sym_tilt()
