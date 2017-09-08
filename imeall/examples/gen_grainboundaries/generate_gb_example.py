from imeall.slabmaker.gengb_from_quat import QuaternionGB
from imeall.slabmaker.slabmaker import build_tilt_sym_gb

#generate a representative array of grain boundary structures:

def gen_canonical_grain_dir(angle, orientation_axis, boundary_plane, material='alphaFe', target_dir='./'):
  """Generate a canonical grain boundary directory. 

  Args:
    angle(float): misorientation angle in radians.
    orienation_axis(list): orientation axis.
    boundary_plane(list): boundary_plane
    material(str, optional): 
    target_dir(str, optional): target directory to deposit canonical grain structure.

  """
  angle_str = str(round((angle*180./np.pi),2)).replace('.', '')
  if len(angle_str) > 4:
    angle_str = angle_str[:-1]
  elif len(angle_str) < 4:
    angle_str = angle_str + '0'

  gbid = ('{0}{1}{2}'.format(orientation_axis[0], orientation_axis[1],orientation_axis[2]) + 
                            angle_str + '{0}{1}{2}'.format(int(abs(gb[1][0])), int(abs(gb[1][1])), int(abs(gb[1][2]))))[0]
  print '\t Grain Boundary ID',  gbid
  target_dir = os.path.join(target_dir, gbid)

  print '\t Grain Boundary Dir', gb_dir
  if not os.path.isdir(target_dir):
    os.mkdir(target_dir)
  else:
    'File already exists.'

  gen_csl(orientation_axis, gb, target_dir=target_dir, gbid=gbid)
  zplanes, dups, nunitcell, grain = build_tilt_sym_gb(gbid, bp=gb[1], v = orientation_axis, 
                                                      c_space=None, target_dir=target_dir,
                                                      rbt=[0.0, 0.0])
  cell = grain.get_cell()
  A    = cell[0][0]*cell[1][1]
  H    = cell[2][2]
  gb_dict = {"gbid"  : gbid, "boundary_plane" : list(gb[1]),
             "orientation_axis" : list(orientation_axis), 
             "type": "symmetric tilt boundary",
             "angle": gb[0], "zplanes" : zplanes, "coincident_sites": dups,
             "n_at" : nunitcell, 'A':A, 'area':A, 'H':H} #area should become dominant keyword

  with open(os.path.join(target_dir, 'gb.json'), 'w') as outfile:
    json.dump(gb_dict, outfile, indent=2)

quat_gb = QuaternionGB()
gb_list = quat_gb.gen_sym_tilt(orientation_axis=[0, 0, 1])

if __name__=='__main__':
#create Canonical Grain Directories from the list.
  for gb in gb_list:
    print gb
    #gen_canonical_grain(gb[0], orientation_axis, gb[1])

