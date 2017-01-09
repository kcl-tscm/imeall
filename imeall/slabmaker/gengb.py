import sys
import numpy as np
import transformations as quat
import operator
from fractions import gcd

orientation_axis = [1., 1., -2.]

#001 axis
if (np.array(orientation_axis) == np.array([0.,0.,1])).all():
  planequat_1      = np.array([0,0,1,0])
  planequat_2      = np.array([0,1,0,0])
#110 axis
elif (np.array(orientation_axis) == np.array([1.,1.,0])).all():
  planequat_1      = np.array([0,0,0,1])
  planequat_1a      = np.array([0,1, -1,2])
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

rad_list = []
deg_list = []

#
#Generate a list of (angle, boundary plane) pairs for
#the chosen orientation axis.
#
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
    rotquat        = quat.quaternion_from_matrix(rotm).round(8)
    rotmquat = (1./(np.array([a for a in rotquat if a!=0.]).min())*rotquat).round(5)
    n1 = quat.quaternion_multiply(planequat_1, rotmquat)
    n2 = quat.quaternion_multiply(planequat_2, rotmquat)
    comm_denom = []
    for a in n1:
      if (a%1).round(1) != 0:
         comm_denom.append(1./(np.abs(a)%1))
    comm_denom = np.array(comm_denom)
    if len(comm_denom)!=0:
     # print '\t {}'.format((comm_denom.max()*n1).round(2))
     # print '\t {}'.format((comm_denom.max()*n2).round(2))
      bp1 = (comm_denom.max()*n1).round(2)
      bp2 = (comm_denom.max()*n2).round(2)
     #bp1 = (reduce(operator.mul, comm_denom)*n1).round(2)
     #bp2 = (comm_denom.max()*n2).round(2)
    else:
      bp1 = n1
      bp2 = n2

    abs_bp1 = map(np.abs, bp1)
    grcd  = reduce(gcd, abs_bp1)
    bp1     = bp1/grcd

    rad  = '\t [{0}, np.array([{1} {2} {3}])],'.format(round(angle, 2),
    bp1[1].round(0), bp1[2].round(0), bp1[3].round(0))
    deg  = '\t [np.pi*({0}/180.), np.array([{1}, {2}, {3}])],'.format(round(180.*angle/np.pi,2), bp1[1].round(2), bp1[2].round(2), bp1[3].round(2))
    deg_list.append([round(180.*angle/np.pi, 2),  [deg]])

deg_list = sorted(deg_list, key=lambda x: x[0])
deg_redun =  [] 
deg_dict  = {}
for deg in deg_list:
  try:
    x=deg_dict[round(deg[0],2)]
  except KeyError:
   deg_dict[round(deg[0],2)]=deg[1][0]
   print deg[1][0]

