from ase.io import read
from ase.io import write
from ase.lattice.spacegroup import crystal
from ase.lattice.surface import surface, bcc111,bcc110
from ase.lattice.cubic   import BodyCenteredCubic
from ase.utils.geometry  import get_duplicate_atoms
from collections import Counter
import json
import os
import numpy as np
import transformations as quat
from quippy import io
from quippy import set_fortran_indexing
try:
  from qlab import set_fortran_indexing, view, gcat
except:
  pass

set_fortran_indexing(False)

def compare_latt_vecs(cell_a, cell_b):
  '''
  Obtain the vector of the ratios of each component of two vectors
  useful for finding lcm of two vectors. Works for orthorhombic cells.
  '''
  ratios = [a[0]/a[1] for a in zip(np.diag(cell_a), np.diag(cell_a))]
  return ratios

def take_pic(fname, translate=False, toggle=False):
  '''
    Rotate and align xyz file for snapshot using AtomEye. 
  '''
  v = view('{0}.xyz'.format(fname))
  v.toggle_bond_mode()
  v.toggle_coordination_coloring()
  v.resize(800,600)
  if translate==True:
    at = gcat()
    z = (at.get_cell()[2,2])/4.
    v.shift_xtal(2, 5.)
  v.rotate([0.,1.,0.],0.5*np.pi)
  v.toggle_parallel_projection()
  # v capture will sometimes toggle parallel projection
  # so double toggle cancels this behaviour.
  if toggle == True:
    v.toggle_parallel_projection()
  v.capture('{0}.png'.format(fname))
  v.close()

def rotate_plane_z(grain, miller):
  ''' 
      Rotates atoms in grain so that planes parallel to plane
      defined by miller index n are parallel to the xy 
      plotting plane.
  '''
  z = np.array([0.,0., 1.])
  p = np.cross(miller, z)
  p *= 1./np.linalg.norm(p)
# Angle between line on plane and z-axis
  theta = np.arccos(miller.dot(z)/(np.linalg.norm(miller)*(np.linalg.norm(z))))
  if theta != theta:
    theta = 0.
    print 'no miller plane specified', theta
    rotation_quaternion = np.array([1., 0., 0., 0.])
  else:
    rotation_quaternion = quat.quaternion_about_axis(theta, p)
  print '\t Aligning Orientation [', int(miller[0]), int(miller[1]), int(miller[2]), '] axis with +z-axis'
  print '\t rotation: ', theta, '(', theta*(180./np.pi), ')'
  print '\t around : ', p, '\n'
  print '\t Rotated Vector (should be along z): ', rotate_vec(rotation_quaternion, miller).round(3)
  grain = rotate_grain(grain, q=rotation_quaternion)
  return rotation_quaternion

def build_tilt_sym_gb(gbid='', bp = [3,3,2], v=[1,1,0],
                      c_space=None, target_dir=None, rbt = None):
  ''' 
      Generate symmetric tilt grain boundary with appropriate configurations: boundary
      plane (bp) oriented along z axis and orthogonal directions in the 
      the other two planes given the orientation axis (v) and an orthogonal vector
      bpxv so we have a proper cube. If rbt is not None then rigid body
      translations are present, this is passed as a list of two numbers
      abs(rbt[0]) < 1. The cell is then displaced as a fxn of these numbers. 
  '''
  bpxv = [(bp[1]*v[2]-v[1]*bp[2]),(bp[2]*v[0]-bp[0]*v[2]),(bp[0]*v[1]- v[0]*bp[1])]
  grain_a = BodyCenteredCubic(directions = [v, bpxv, bp],
                              size = (1,1,1), symbol='Fe', pbc=(1,1,1),
                              latticeconstant = 2.83)

  n_grain_unit = len(grain_a)
# For the Chamati EAM potential the cutoff radius is 5.67 A
# We want to separate the grain boundaries by at least this much, A safe
# minimum value is 12 A. 
  n = 2
  while(grain_a.get_cell()[2,2]< 12.0 and n < 10):
    grain_a = BodyCenteredCubic(directions = [v, bpxv, bp],
                              size = (1,1,n), symbol='Fe', pbc=(1,1,1),
                              latticeconstant = 2.83)
    n += 1
  print '\t {0} repeats in z direction'.format(n)
  grain_b = grain_a.copy()
  grain_c = grain_a.copy()
  print '\t', '{0} {1} {2}'.format(v[0],v[1],v[2])
  print '\t', '{0} {1} {2}'.format(bpxv[0],bpxv[1],bpxv[2])
  print '\t', '{0} {1} {2}'.format(bp[0], bp[1], bp[2])
  if c_space==None:
    s1 = surface('Fe', (map(int, bp)), n)
    c_space = s1.get_cell()[2,2]/float(n) #-s1.positions[:,2].max()
    s2 = surface('Fe', (map(int, v)), 1)
    x_space = s2.get_cell()[0,0] #-s1.positions[:,2].max()
    s3 = surface('Fe', (map(int, bpxv)), 1)
    y_space = s3.get_cell()[1,1] #-s1.positions[:,2].max()
  print '\t Interplanar spacing: ', x_space.round(2), y_space.round(2), c_space.round(2), 'A'
# Reflect grain b in z-axis (across mirror plane):
  print grain_a.get_cell()[2,2]-grain_a.positions[:,2].max()
  grain_b.positions[:,2]  = -1.0*grain_b.positions[:,2]
  grain_c.extend(grain_b)
  grain_c.set_cell([grain_c.get_cell()[0,0], grain_c.get_cell()[1,1], 2.*grain_c.get_cell()[2,2]])
  grain_c.positions[:,2] += abs(grain_c.positions[:,2].min())
  pos = [grain.position for grain in grain_c]
  pos = sorted(pos, key= lambda x: x[2])
  dups = get_duplicate_atoms(grain_c)
# Now build the second grain boundary by reflecting in y plane
  grain_b = grain_c.copy()
  grain_b.positions[:,0] = -grain_b.positions[:,0]
  grain_b.positions[:,1] = -grain_b.positions[:,1]
  grain_b.wrap()
# In BCC lattice spacing is different depending on whether the
# z-plane has an even or odd number of odd miller indices:
  if sum([int(n)%2 for n in bp])%2 == 0 :
    grain_b.positions[:,2] -= (grain_b.positions[:,2].max()+2*c_space)
  else:
    grain_b.positions[:,2] -= (grain_b.positions[:,2].max()+c_space)
  grain_c.extend(grain_b)
  grain_c.set_cell([grain_c.get_cell()[0,0], grain_c.get_cell()[1,1], 2.*grain_c.get_cell()[2,2]])
  grain_c.positions[:,2] += abs(grain_c.positions[:,2].min())
  dups = get_duplicate_atoms(grain_c)
# Displace replicated atoms along the unit cell
  for dup in dups:
      grain_c[dup[1]].position[2] += grain_a.get_cell()[2,2]
  shared_atoms = [grain_c[dup[0]] for dup in dups]
  for at in shared_atoms:
    print '\t', at.index, at.position
  print '\t There are {0} shared atoms'.format(len(shared_atoms))
  z_planes = [round(atom.position[2],4) for atom in shared_atoms]
  z_planes = list(sorted(set(z_planes)))
  cell = grain_c.get_cell()
  if rbt != None:
    for grain in grain_c:
      if (z_planes[0] < grain.position[2].round(4) < z_planes[1]):
        grain.position[0] += rbt[0]*cell[0][0]
        grain.position[1] += rbt[1]*cell[1][1]
#so move one of the coincident atoms if there are any.
    for atom in shared_atoms:
      grain_c.positions[atom.index, 0] += rbt[0]*cell[0][0]
      grain_c.positions[atom.index, 1] += rbt[1]*cell[1][1]
    grain_c.wrap()
#Grain was explicitly set to 2*grain A before so we need 
#to space out an additional layer.
  if sum([int(n)%2 for n in bp])%2 == 0 :
	  grain_c.center(vacuum=c_space/2., axis=2)
  else:
    grain_c.center(vacuum=c_space/4., axis=2)
# if we specify a target directory take a picture and
# write the file there directly, otherwise just return the grain
# for further processing.
  if target_dir != None:
    print '\t Writing {0}.xyz to file'.format(gbid)
    io.write('{0}.xyz'.format(os.path.join(target_dir, gbid)), grain_c)
    try:
      take_pic(os.path.join(target_dir, gbid),toggle=True)
    except:
      'atomeye not pleased!'
    return [z_planes, len(dups), n_grain_unit, grain_c]
  else:
    return grain_c

def rotate_plane_y(grain, miller):
  # rotates atoms in grain so that the miller plane 
  # in grain is parallel to the y axis.
  # Returns the quaternion characterizing the rotation which accomplishes this.
  y = np.array([0., 1., 0.])
  p = np.cross(miller, y)
  p *= 1./np.linalg.norm(p)
  # angle between line in plane and desired line:
  theta = np.arccos(miller.dot(y)/(np.linalg.norm(miller)*(np.linalg.norm(y))))
  if theta != theta:
    theta = 0.
    rotation_quaternion = np.array([1., 0., 0., 0.])
    print 'no miller plane specified', theta
  else:
    rotation_quaternion = quat.quaternion_about_axis(theta, p)
    print '\t quaternion: ', rotation_quaternion.round(2)
  print '\t Rotated miller index: ', miller.round(3)
  print '\t Aligning miller plane with y-axis by angle: ', theta*(180./np.pi)
  print '\t around : ', p.round(3)
  print '\t Rotated vector (should be along y): ', rotate_vec(rotation_quaternion, miller).round(3), '\n'
  grain = rotate_grain(grain, q = rotation_quaternion)
  return rotation_quaternion

def generate_mirror(grain, point, miller):
#reflection matrix for plane around point
#with miller index as normal to that plane
  M = quat.reflection_matrix(point, miller) 
  print 'mirror plane', M
  for atom in grain:
    p = atom.position
    p_prime = M[:3,:3].dot(p)
    atom.position = p_prime
#
#  From  zeiner we have the following theorem and wisdom:
#  Firstly a lattice is coincident if and only if it 
#  is an orthogonal matrix with rational entries.
#  The general relationship between a rotation axis, rotation, and a 
#  quaternion can be written:
#  In order to be pure we should !!not!! do the rotations in terms of !!Degrees!!
#  It is superior to do the rotation in terms of an integral number of radians.
#  Then we can exploit all of Zeiner's relationships. For instance:
#  quaternion rotation angle \psi and the lattice 
#  direction of the rotation axis.
#  Additionally the coincidence index can be given by
#  m is angle of rotation in radians, atoms is the grain to be rotated
# (m,n,0,0) := \psi = \arccos(\frac{m^{2} - n^{2}}{m^{2} + n^{2}}), [100]
#
def find_sigma_csl(q):
# \Sigma(R(\mathbf{r})) = |\mathbf{r}|^{2}/2^{l} where l is the maximal power such that
# 2^{l} divides |\mathbf{r}|^{2}.
# m is angle of rotation in radians
#(m,n,n,0) := \psi = \arccos(\frac{m^{2} - 3n^{2}}{m^{2} + 3n^{2}}), [111]
#(k,lambda,mu,nu)   :=
#       \psi = \arccos(\frac{k^{2} - \lambda^{2} - \mu^{2} -\nu^{2}}{k^{2} +
#           +  \lambda^{2} + \mu^{2} + \nu^{2}}), [\lambda\mu\nu]
  r   = q.dot(q)
  div = r
  l   = 0
  while div >= 1:
    div = r/np.power(2,l)
    l += 1
  return sigma

def gnu_plot_gb(boundary_plane, m, invm, gbid, mb=0.0, invmb=0.0,target_dir='./'):
# output gnuplot script to generate "CSL plots on the fly"
# in svg format for visualization in imeall.
  f = open(os.path.join(target_dir, 'plot.gnu'), 'w')
  script = " \
  h(x) = {3}*x   \n \
  set xr[-10:10] \n \
  set yr[-15:15] \n \
  pl   'grainaT.dat' u 1:2 w p pt 7 ps 1.3 lt 1 t 'Grain A' \n \
  repl 'grainaB.dat' u 1:2 w p pt 6 lt 1.0      t ''        \n \
  repl 'grainbT.dat' u 1:2 w p pt 7 ps 1 lt 3 t 'Grain B'   \n \
  repl 'grainbB.dat' u 1:2 w p pt 6 ps 1 lt 3 t ''          \n \
  repl h(x) lt -1 t '({0})' \n \
  set terminal svg          \n \
  set output 'csl_{1}.svg'  \n \
  repl                      ".format(boundary_plane, gbid, m, invm, mb, invmb)
  print >> f, script
  f.close()

def bcc_csl_nn0(m, n, grain):
# m
#(m,n,n,0) := \psi = \arccos(\frac{m^{2} - 2n^{2}}{m^{2} + 2n^{2}}), [110]
# once m and n have been chosen according to some scheme we generate the rotation
# matrix
# R = quat.quaternion_matrix([m,n,n,0])
  R = zeiner_matrix(np.array([m,n,n,0]))
# project down to 3x3:
  for i in range(len(grain)):
    grain[i].position = np.dot(R, grain[i].position)

def zeiner_matrix(q):
  r2 = q.dot(q)
  k = q[0]
  l = q[1]
  m = q[2]
  n = q[3]
  R = np.array([
        [ k*k+l*l-m*m-n*n,   -2*k*n + 2*l*m,    2*k*m + 2*l*n],
        [   2*k*n + 2*l*m,  k*k-l*l+m*m-n*n, -2*k*l + 2*m*n],
        [  -2*k*m + 2*l*n, 2*k*l + 2*m*n,  k*k-l*l-m*m + n*n]
        ])
# The zeiner matrix refers to the rotation matrix
# for the lattice when applied to grains we require:
# R.T
  R = (1./r2)*R.T
  return R

def csl_lattice_vecs(m,n):
# rational orthogonal matrices can be parameterized
# by integral quaternions, i.e. by quaternions with integral coefficients or
# a quaternion with integral coefficients + 1/2(1 1 1 1).
# an integral quaternion is primitive if the greatest divisor of its integral components is 1.
# for a primitive quaternion \mathbf{r} = (r_0,r_1,r_2,r_3)
# r^{0} = [ r_{1}, r_{2}, r_{3}]
# r^{1} = [r_{0}, r_{3}, -r_{2}]
# r^{2} = [-r_{3}, r_{0}, r_{1}]
# r^{3} = [r_{2}, -r_{1}, r_{0}]
# Given a bcc lattice (Zeiner05) we can define the
# lattice vectors of the coincident site lattice as
# r^{0}, r^{1}, r^{2}, r^{3}, 1/2(r^{0} + r^{1} + r^{2} + r^{3}),  if |r|^{2} is odd
# r^{0}, 1/2(r^{0}+r^{1}), 1/2(r^{0} + r^{2}), 1/2(r^{0} + r^{3}), if 2 divides |r|^2
# and 4 does not divide |r|^{2}
# 1/2 r^{0}, 1/2 r^{1}, 1/2r^{2}, 1/2 r^{3} if 4 divides |r|^2
# There is then a straight forward procedure to obtain a 2d lattice
# for the grain boundary. Just take the miller plane defining the boundary
# and pick a set of independent vectors from the lattice vectors of the CSL.
  r = np.array([m,n,n,0])
  r0 = np.array([r[1], r[2],  r[3]])
  r1 = np.array([r[0], r[3], -r[2]])
  r2 = np.array([-r[3], r[0], r[1]])
  r3 = np.array([r[2], -r[1], r[0]])
  if np.mod(r.dot(r),2)==1:
    print '\t |r|^{2} is odd'
    a = r0
    b = r1
    c = r2
    d = 0.5*(r0 + r1+ r2 +r3)
  elif (np.mod(r.dot(r), 2) == 0 and np.mod(r.dot(r),4) != 0):
    print '\t |r|^{2}  is divisible by 2 but not 4'
    a = r0
    b = 0.5*(r0 +r1)
    c = 0.5*(r0+r2)
    d = 0.5*(r0 + r3)
# elif mod(r.dot(r), 4)==0:
  else:
    print '\t |r|^{2}  is divisible 4'
    a = 0.5*r0
    b = 0.5*(r1)
    c = 0.5*(r2)
    d = 0.5*(r3)
  return np.array([a,b,c,d])

def rotate_vec(q, vec):
  '''
    Rotate 3 vector with quaternion
  '''
  vec = np.array([0., vec[0], vec[1], vec[2]])
  qm = quat.quaternion_conjugate(q)
  pos_prime = quat.quaternion_multiply(q, vec)
  pos_prime = quat.quaternion_multiply(pos_prime, qm)
  vec = quat.quaternion_imag(pos_prime)
  return vec

def rotate_grain(grain, theta=0., x=[0.,0.,0.], q=[]):
  '''
    Rotate the grain according to quaternion = [Theta, u,v,w] = [w,x,y,z]
    Standard routine is passed angle and vector
    this generates the quaternion to do the rotation 
    however if a quaternion, q!=None, is passed to routine the plane is rotated
    using q.
  '''
  if q==[]:
    q = quat.quaternion_about_axis(theta, x)
    qm = quat.quaternion_conjugate(q)
  else:
    qm = quat.quaternion_conjugate(q)
  for i in range(len(grain)):
    pos       = np.array([0, grain[i].position[0], grain[i].position[1], grain[i].position[2]])
    pos_prime = quat.quaternion_multiply(q, pos)
    pos_prime = quat.quaternion_multiply(pos_prime, qm)
    grain[i].position = quat.quaternion_imag(pos_prime)
  return grain

def print_points(atoms, f):
  for atom in atoms:
    print>>f, '{0:12.7f} {1:12.7f} {2:12.7f}'.format(atom.position[0], atom.position[1], atom.position[2])

def find_densest_plane(grain_dict):
  maxx = max([len(a) for a in grain_dict.values()])
  len_keys = {x:len(y) for x,y in grain_dict.items()}
  keys_2   = [x for x,y in grain_dict.items()]
  keys_2   = sorted(keys_2)
  keys = [x for x,y in grain_dict.items() if len(y) == maxx]
  print '\t Largest number of points: ', maxx, 'Number of z-planes: ', len(keys)
  key1 = keys[0]
  z_below = keys_2.index(key1)-1
  key2 = keys_2[z_below]
  return key1, key2

def simplify_csl(m, b=0.00,target_dir='./'):
# Simplify the dat files of the csl so that atoms in grain_a are below the line
# defined by the boundary plane in the current projection and atoms from grain_b
# are above.
  graina_list = ['grainaT.dat', 'grainaB.dat']
  grainb_list = ['grainbT.dat', 'grainbB.dat']
  for name in graina_list:
    f = open('{0}'.format(os.path.join(target_dir, name)), 'r+')
    grain_a = f.read().split('\n')
    f.close()
    f = open('{0}'.format(os.path.join(target_dir, name)), 'w')
    del(grain_a[-1])
    for at in grain_a:
      position = np.array(map(float, at.split()))
      if ( round(position[1]- m*position[0]-b,2) <= 0.0):
        print>>f, '{0:12.7f} {1:12.7f} {2:12.7f}'.format(position[0], position[1], position[2])
      else:
        pass
    f.close()
  for name in grainb_list:
    f = open('{0}'.format(os.path.join(target_dir, name)),'r+')
    grain_b = f.read().split('\n')
    f.close()
    f = open('{0}'.format(os.path.join(target_dir, name)),'w')
    del(grain_b[-1])
    for at in grain_b:
      position = np.array(map(float, at.split()))
      if (round(position[1]-m*position[0]-b,2) >= 0.0):
        print>>f, '{0:12.7f} {1:12.7f} {2:12.7f}'.format(position[0], position[1], position[2])
      else:
        pass
    f.close()

def gen_csl(orientation_axis, gb, target_dir='./', gbid='0000000000'):
    a  = 2.834
    fe = crystal('Fe', [(0,0,0)], spacegroup=229,
                 cellpar=[a, a, a, 90, 90, 90],size=[10,10,10])
    grain_a = fe.copy()
    grain_b = fe.copy()
    for i in range(3):
      grain_a.positions[:,i] -= 10*a/2.
      grain_b.positions[:,i] -= 10*a/2.
    if np.allclose(orientation_axis, [1,1,0]):
      n = gb[1][0]
      m = gb[1][2]
    elif np.allclose(orientation_axis, [0,0,1]):
      n = gb[1][0]
      m = gb[1][1]
    elif np.allclose(orientation_axis, [1,1,1]):
      m = gb[1][0]
      n = gb[1][1]
    rotm           = quat.rotation_matrix(gb[0], orientation_axis)
    rotquat        = quat.quaternion_from_matrix(rotm).round(8)
    angle_str      = str(round((gb[0]*180./np.pi),2)).replace('.', '')
#   Truncate the angle string if it is greater than 4 characters
    if len(angle_str) > 4:
      angle_str = angle_str[:-1]
    elif len(angle_str) < 4:
      angle_str = angle_str + '0'
    print '\t Grain Boundary ID: ', gbid
    print '\t Rotation Axis: ', orientation_axis, 'Rotation Angle: ', (round(gb[0],6))*(180./np.pi),round(gb[0],4)
    if (orientation_axis == [1,1,0]).all():
      print '\t Angle from Zeiner Formula', np.arccos((m*m-2.*n*n)/(m*m+2.*n*n))*(180./np.pi), m*m-2*n*n, m*m+2*n*n
    elif (orientation_axis == [0,0,1]).all():
      print '\t Orientation axis  0 0 1 '
      print '\t Angle from Zeiner Formula', np.arccos((m*m-n*n)/(m*m+n*n))*(180./np.pi), m*m-n*n, m*m+n*n
    elif np.allclose(orientation_axis, [1,1,1]):
      print '\t Orientation axis: ', orientation_axis
      print '\t Angle from Zeiner Formula', np.arccos((m*m-3*n*n)/(m*m+3*n*n))*(180./np.pi), m*m-3*n*n, m*m+3*n*n

    print '\t Rotation quaternion: ', rotquat
    print '\t Integral Rotation quaternion: ', (1./(np.array([a for a in rotquat if a !=0]).min())*rotquat).round(5)
    print '\t Boundary plane grain A coordinate system: ', gb[1].round(3)

    rotm           = quat.rotation_matrix(gb[0], orientation_axis)
    rotmquat       = (1./(np.array([a for a in rotquat if a !=0.]).min())*rotquat).round(5)

    if np.allclose(orientation_axis, [1,1,0]):
      planequat_1    = np.array([0,0,0,1])  
      planequat_2    = np.array([0,1,-1,0]) 
    elif np.allclose(orientation_axis, [0,0,1]):
      planequat_1    = np.array([0,0,1,0])  
      planequat_2    = np.array([0,1,0,0]) 
    elif np.allclose(orientation_axis, [1,1,1]):
      planequat_1    = np.array([0,1,-1,0])  
      planequat_2    = np.array([0,1,1,-2]) 

    print '\t These vectors satisfy the conditions for a Symmetric Boundary Plane: '
    n1 = quat.quaternion_multiply(planequat_1, rotmquat)
    n2 = quat.quaternion_multiply(planequat_2, rotmquat)
    comm_denom = []
    for a in n1:
      if (a%1).round(1) != 0:
        comm_denom.append(1./(np.abs(a)%1))
    comm_denom = np.array(comm_denom)
    if len(comm_denom)!=0:
      print '\t {}'.format((comm_denom.max()*n1).round(2))
      print '\t {}'.format((comm_denom.max()*n2).round(2))
      print '\n'
      print '\n'
    else:
      print '\t', n1
      print '\t', n2
      print '\n'
      print '\n'
    theta          = gb[0]
    boundary_plane = gb[1]
    if theta != 0:
      m2 = 2*(1.+np.cos(theta))/(1.-np.cos(theta))
    else:
      m2 = 0
    if theta !=0:
      m = np.sqrt(m2).round(0)
      n = 1
    else:
      m=0
      n=0

    csl_factory(orientation_axis, boundary_plane, m, n, grain_a, grain_b,
                theta=theta, mode='', target_dir= target_dir)

def csl_factory(orientation_axis, boundary_plane, m, n, grain_a, grain_b,
                mode='Zeiner', theta=0, target_dir='./'):
  """ 
      Builds a plot of the coincident site lattice oriented along a suitable
      plane for viewing purposes. In the gnuplot plane we think of x,y,z
      General orientation axis = N = z
      Vector specifying 0 Angle (for the case of tilt) = v = x
      y = Nxv. rotate_plane_z takes grain_a and rotates it
      so that the orientation axis of the grain boundary with
      respect to grain a is orthogonal to the x-y plane. The
      quaternion to accomplish this rotation is stored in plane_quaternion_z 
  """
  f = open(os.path.join(target_dir, 'grainaT.dat'), 'w')
  g = open(os.path.join(target_dir, 'grainaB.dat'), 'w')
  if np.allclose(orientation_axis, [0,0,1]):
    print 'Axis already aligned along z'
    plane_quaternion_z = np.array([1,0,0,0])
  else:
    plane_quaternion_z = rotate_plane_z(grain_a, (-1.*orientation_axis))

  print '\t boundary plane', boundary_plane.round(2)
# plane_quaternion_x stores the rotation quaternion so that the boundary plane in grain_a coordinates
# can be rotated so that it is perpendicular to the x-axis.
  plane_quaternion_x = rotate_plane_y(grain_a, rotate_vec(plane_quaternion_z, boundary_plane))
  print '\t boundary plane', boundary_plane.round(2)
# The following code finds the planes with largest number
# of points for viewing purposes:
# First we sort grain_a by x, y, z
  grain_a = sorted(grain_a, key=lambda x: (x.position[2], x.position[0], x.position[1]))
  z_vals = []
# The loop over each position in grain_a and store the z-value
# in a list.
  for grain in grain_a:
    if round(grain.position[2], 2) not in [round(z,2) for z in z_vals[:]]:
      z_vals.append(round(grain.position[2],2))
  grain_dict = {}
# initialize dictionary of lists
  for z in z_vals:
    grain_dict[z] = []
# append the atoms to each position in the dictionary.
  for grain in grain_a:
    grain_dict[round(grain.position[2], 2)].append(grain)
# create a list of z-values and the number of points for each z value
  len_list = []
  for z in z_vals:
     len_list.append((z, len(grain_dict[z])))
# keys = [x for x,y in grain_dict.items() if len(y) == maxx]
  key1, key2 = find_densest_plane(grain_dict)
  print_points(grain_dict[key1], f)
  print_points(grain_dict[key2], g)
  f.close()
  g.close()
# Setup grain b and project into chosen plane
  f = open(os.path.join(target_dir, 'grainbT.dat'),'w')
  g = open(os.path.join(target_dir, 'grainbB.dat'),'w')
  if mode=='Zeiner':
    print '\t Using Zeiner'
    print '\t Spanning Vectors of CSL'
    print '\t', csl_lattice_vecs(m,n)
    print ''
    if (m!=0 and n!=0):
      bcc_csl_nn0(m, n, grain_b)
# A rotation is a coincidence rotation if the rotation
# matrix is orthogonal and has rational entries.
    if (m!=0 and n !=0):
      rotm = zeiner_matrix(np.array([m,n,n,0]))
      print rotm
    else:
      rotm = np.identity(3)
    if m!=0 and n!=0:
      mn = float(m*m-2*n*n)/float(m*m+2*n*n)
    else:
      mn = 0
    angle_of_rotation = np.arccos(mn)
    deg = round(angle_of_rotation*(180./np.pi),2) 
    print m, n
    print '\n Angle of rotation: ', deg, 'or ', 180.-deg, ' Orientation Axis: \n', (n, n, 0), '\n'
  else:
# Theta = np.arccos(float(m-2*n*n)/float(m+2*n*n))
    rotate_grain(grain_b, -theta, orientation_axis)
    deg = round(theta*(180./np.pi), 2) 
    print '\n Angle of rotation: ', deg, 'or ', 180.-deg, ' Orientation Axis: ', (n, n, 0), '\n'

  intercept = np.array([-1.0, 1.0, 2.0])
  boundary_plane = rotate_vec(plane_quaternion_z, boundary_plane)
  boundary_plane = rotate_vec(plane_quaternion_x, boundary_plane)
  intercept = rotate_vec(plane_quaternion_z, intercept)
  intercept = rotate_vec(plane_quaternion_x, intercept)
  print ''
  print 'Boundary plane: ', boundary_plane.round(3)
  print 'Intercept: ', intercept.round(3)
  print 'Slope of line in x,y plane'
# Draw a line on the graph which is the 
# equation of the line defining the miller index
# and a line perpendicular to this which corresponds
# to a line in the plane.
# m = (boundary_plane[0]/boundary_plane[1])
# invm = -(1./(boundary_plane[0]/boundary_plane[1]))
# Grain boundary should be parallel to \hat{x}
  m    = 0.
  invm = 0.
  print '   m = ', m
  print 'invm = ', invm
#
# The boundary_plane vector is the miller index
# of the grain boundary with respect to initial crystal
# for grain b we rotate this around orientation axis
# by the orientation angle and then.
#
  boundary_plane_b = boundary_plane
  rotm = zeiner_matrix(np.array([m,n,n,0]))
  boundary_plane_b = np.dot(rotm, boundary_plane_b)
  print '\t Rotated Grain Boundary:' 
  boundary_plane_b = rotate_vec(plane_quaternion_z, boundary_plane_b)
  boundary_plane_b = rotate_vec(plane_quaternion_x, boundary_plane_b)
  print '\t', boundary_plane_b
  mb    = 0.
  invmb = 0.
  gnu_plot_gb('Plane', m, invm, gbid, mb=mb, invmb=invmb, target_dir=target_dir)
  rotate_grain(grain_b, q=plane_quaternion_z)
  rotate_grain(grain_b, q=plane_quaternion_x)
  grain_b = sorted(grain_b, key=lambda x: (x.position[2], x.position[0], x.position[1]))
  z_vals = []
  for grain in grain_b:
    if round(grain.position[2], 2) not in [round(z,2) for z in z_vals[:]]:
      z_vals.append(round(grain.position[2],2))
  grain_dict = {}
  for z in z_vals:
    grain_dict[z] = []
  for grain in grain_b:
    grain_dict[round(grain.position[2], 2)].append(grain)
  len_list = []
  for z in z_vals:
     len_list.append((z,len(grain_dict[z])))
# maxx = max([len(a) for a in grain_dict.values()]) 
# keys = [x for x,y in grain_dict.items() if len(y) == maxx ]
  key1, key2 = find_densest_plane(grain_dict)
  print_points(grain_dict[key1], f)
  print_points(grain_dict[key2], g)
  f.close()
  g.close()
  simplify_csl(0.,target_dir=target_dir)

# BCC Iron unit cell
if __name__=='__main__':
# These lists shouldn't be necessary if we exploit
# the relationship between the rotation quaternion
# and the tilt symmetric mirror plane.
# Angle from zeiner
# 100, \psi = ( m**2 -   n**2)/(m**2-n**2)
# 110, \psi = ( m**2 - 2*n**2)/(m**2-n**2)
# 111, \psi = ( m**2 - 3*n**2)/(m**2-n**2)
  sym_tilt_100 = [[np.pi*(53.13/180.),  np.array([2,1,0])],   
                  [np.pi*(36.87/180.),  np.array([3,1,0])],
                  [np.pi*(28.07/180.),  np.array([4,1,0])],
                  [np.pi*(22.62/180.),  np.array([5,1,0])],
                  [np.pi*(18.92/180.),  np.array([6,1,0])],
                  [np.pi*(16.26/180.),  np.array([7,1,0])],
                  [np.pi*(14.25/180.),  np.array([8,1,0])],
                  [np.pi*(12.68/180.),  np.array([9,1,0])],
                  [np.pi*(11.42/180.),  np.array([10,1,0])],
                  [np.pi*(67.38/180.),  np.array([3,2,0])],
                  [np.pi*(53.13/180.),  np.array([4,2,0])],
                  [np.pi*(43.60/180.),  np.array([5,2,0])],
                  [np.pi*(36.87/180.),  np.array([6,2,0])],
                  [np.pi*(31.89/180.),  np.array([7,2,0])],
                  [np.pi*(28.07/180.),  np.array([8,2,0])],
                  [np.pi*(25.05/180.),  np.array([9,2,0])],
                  [np.pi*(22.61/180.),  np.array([10,2,0])],
                  [np.pi*(143.13/180.), np.array([1,3,0])],   
                  [np.pi*(73.739/180.), np.array([4,3,0])],   
                  [np.pi*(61.93/180.),  np.array([5,3,0])],   
                  [np.pi*(53.13/180.),  np.array([6,3,0])],   
                  [np.pi*(46.40/180.),  np.array([7,3,0])],   
                  [np.pi*(41.11/180.),  np.array([8,3,0])],   
                  [np.pi*(36.87/180.),  np.array([9,3,0])],   
                  [np.pi*(33.39/180.),  np.array([10,3,0])],   
                  [np.pi*(81.20/180.),  np.array([7,6,0])],
                  [np.pi*(104.25/180.), np.array([7,9,0])],
                  [np.pi*(119.48/180.), np.array([7,12,0])],
                  [np.pi*(47.92/180.),  np.array([4, 9, 0])]
                  ]
# 0 degrees is at [1-10] rotation by 180 would set
# the zero of the angle as [-110].
  sym_tilt_110 = [[np.pi*(13.44/180.), np.array([-1., 1., 12.])],
                  [np.pi*(20.05/180.), np.array([-1., 1., 8.])],
                  [np.pi*(26.53/180.), np.array([-1., 1., 6.])],
                  [np.pi*(31.59/180.), np.array([-1., 1., 5.])],
                  [np.pi*(38.94/180.), np.array([-1., 1., 4.])],
                  [np.pi*(44.00/180.), np.array([-2., 2., 7.])],
                  [np.pi*(50.48/180.), np.array([-1., 1., 3.])],
                  [np.pi*(58.99/180.), np.array([-2., 2., 5.])],
                  [np.pi*(70.53/180.), np.array([-1, 1, 2])],
                  [np.pi*(80.63/180.), np.array([-3, 3, 5])],
                  [np.pi*(86.63/180.), np.array([-2, 2, 3])],
                  [np.pi*(93.37/180.), np.array([-3., 3., 4.])],
                  [np.pi*(99.37/180.), np.array([-3., 3., 4.])],
                  [np.pi*(109.47/180.), np.array([-1., 1., 1.])],
                  [np.pi*(121.01/180.), np.array([-5., 5., 4.])],
                  [np.pi*(129.52/180.), np.array([-3., 3., 2.])],
                  [np.pi*(141.06/180.), np.array([-2., 2., 1.])],
                  [np.pi*(148.41/180.), np.array([-5., 5., 2.])],
                  [np.pi*(153.47/180.), np.array([-3., 3., 1.])],
                  [np.pi*(159.95/180.), np.array([-4., 4., 1.])],
                  [np.pi*(166.55/180.), np.array([-6., 6., 1.])]
                 ]

  sym_tilt_111 = [ 
				   [np.pi*(19.65/180.), np.array([9.0, -11.0, 2.0])],
				   [np.pi*(19.65/180.), np.array([-13.0, -7.0, 20.0])],
				   [np.pi*(21.79/180.), np.array([4.0, -5.0, 1.0])],
				   [np.pi*(21.79/180.), np.array([-2.0, -1.0, 3.0])],
				   [np.pi*(24.43/180.), np.array([7.0, -9.0, 2.0])],
				   [np.pi*(24.43/180.), np.array([-11.0, -5.0, 16.0])],
				   [np.pi*(27.8/180.),  np.array([3.0, -4.0, 1.0])],
				   [np.pi*(27.8/180.),  np.array([-5.0, -2.0, 7.0])],
				   [np.pi*(32.2/180.),  np.array([5.0, -7.0, 2.0])],
				   [np.pi*(32.2/180.),  np.array([-3.0, -1.0, 4.0])],
				   [np.pi*(38.21/180.), np.array([2.0, -3.0, 1.0])],
				   [np.pi*(38.21/180.), np.array([-4.0, -1.0, 5.0])],
				   [np.pi*(42.1/180.),  np.array([2.0, -3.0, 1.0])],
				   [np.pi*(46.83/180.), np.array([3.0, -5.0, 2.0])],
				   [np.pi*(46.83/180.), np.array([-7.0, -1.0, 8.0])],
				   [np.pi*(60.0/180.),  np.array([1.0, -2.0, 1.0])],
				   [np.pi*(75.18/180.), np.array([1.0, -3.0, 2.0])],
				   [np.pi*(81.79/180.), np.array([1.0, -3.0, 2.0])],
				   [np.pi*(81.79/180.), np.array([-5.0, 1.0, 4.0])],
				   [np.pi*(89.41/180.), np.array([1.0, -3.0, 2.0])],
				   [np.pi*(89.41/180.), np.array([-5.0, 1.0, 4.0])],
				   [np.pi*(92.2/180.),  np.array([1.0, -3.0, 2.0])],
				   [np.pi*(108.36/180.), np.array([-2.0, 1.0, 1.0])],
				   [np.pi*(120.0/180.),  np.array([-2.0, 1.0, 1.0])],
				   [np.pi*(133.17/180.), np.array([-2.0, 1.0, 1.0])],
				   [np.pi*(137.9/180.), np.array([-3.0, 2.0, 1.0])],
				   [np.pi*(147.8/180.), np.array([-1.0, -3.0, 4.0])],
				   [np.pi*(147.8/180.), np.array([-7.0, 5.0, 2.0])],
				   [np.pi*(158.21/180.), np.array([-1.0, -2.0, 3.0])],
				   [np.pi*(158.21/180.), np.array([-5.0, 4.0, 1.0])],
				   [np.pi*(163.57/180.), np.array([-3.0, -5.0, 8.0])],
           [np.pi*(163.57/180.), np.array([-13.0, 11.0, 2.0])]
              ]
#Orphan Boundaries:
  test = [[np.pi*(92.2/180.),  np.array([2.0, 3.0, -5.0])]]

#These do not appear to form symmetric tilt boundaries
  Random = [[np.pi*(92.2/180.),  np.array([-5.0, 1.0, 3.0])],
				  [np.pi*(104.82/180.), np.array([-4.0, 2.0, 3.0])]]
#for gb in sym_tilt_110[:]:
#orientation_axis = np.array([0, 0, 1])
#orientation_axis = np.array([1, 1, 0])
#for gb in sym_tilt_110[8:10]:
  orientation_axis = np.array([1, 1, 1])
  for gb in test:
  #for gb in sym_tilt_111[:]:
    angle_str      = str(round((gb[0]*180./np.pi),2)).replace('.', '')
    if len(angle_str) > 4:
      angle_str = angle_str[:-1]
    elif len(angle_str) < 4:
      angle_str = angle_str + '0'
    gbid = '{0}{1}{2}'.format(orientation_axis[0], orientation_axis[1], 
        orientation_axis[2]) + angle_str + '{0}{1}{2}'.format(int(abs(gb[1][0])), int(abs(gb[1][1])), int(abs(gb[1][2])))
    print '\t Grain Boundary ID',  gbid
#Dump GBs in this directory:
    gb_dir     = os.path.join('./alphaFeDFT','111')
    target_dir = os.path.join(gb_dir, gbid)
    print '\t Grain Boundary Dir', gb_dir
    if not os.path.isdir(target_dir):
      os.mkdir(target_dir)
    else:
      'directory already exists'
# Drop the coincident site lattice stuff in the directory
    gen_csl(orientation_axis, gb, target_dir=target_dir, gbid=gbid)
# Drop the grain.xyz file in the target directory
    zplanes, dups, nunitcell, grain_c = build_tilt_sym_gb(gbid, bp=gb[1], v = orientation_axis, 
                                                 c_space=None,
                                                 target_dir=target_dir,
                                                 rbt=[0.0, 0.75])
# json id file contains, gbid, boundaryplane, zplanes(coordinates of center of
# grain boundary
    cell = grain_c.get_cell()
    A    = cell[0][0]*cell[1][1]
    H    = cell[2][2]
    gb_dict = { "gbid"  : gbid, "boundary_plane" : list(gb[1]),
                "orientation_axis" : list(orientation_axis), 
                "type": "symmetric tilt boundary",
                "angle": gb[0], "zplanes" : zplanes, "coincident_sites": dups,
                "n_at" : nunitcell, 'A':A , 'H':H}

    with open(os.path.join(target_dir, 'gb.json'), 'w') as outfile:
      json.dump(gb_dict, outfile, indent=2)
    ovito = "~/ovito-2.6.1-x86_64/bin/ovito"
    if(not True):
      if os.path.isfile(ovito):
        os.system("{0} {1}.xyz".format(ovito, os.path.join(target_dir, gbid)))
        variable = raw_input('Continue?')
      elif os.path.isfile('/Users/lambert/Ovito.app/Contents/MacOS/ovito'):
        os.system('/Users/lambert/Ovito.app/Contents/MacOS/ovito {0}.xyz'.format(os.path.join(target_dir, gbid)))
        variable = raw_input('Continue?')
      else:
        print 'ovito is in non-standard place cannot display crystal structure'
        variable = raw_input('Continue?')
      if variable == 'y':
        pass
      elif variable =='n':
        break
