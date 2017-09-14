import os
import sys
import json
import numpy as np
from imeall.slabmaker.gengb_from_quat import QuaternionGB
from imeall.slabmaker.slabmaker import build_tilt_sym_gb, gen_csl, gen_canonical_grain_dir

#generate a representative array of grain boundary structures:

quat_gb = QuaternionGB()
or_axis = [0,0,1]
or_axis = [1,1,0]
#or_axis = [1,1,1]
gb_list = quat_gb.gen_sym_tilt(orientation_axis=or_axis)

if __name__=='__main__':
#create Canonical Grain Directories from the list.
  for gb in gb_list[:5]:
    print gb[0], gb[1]
    gen_canonical_grain_dir(gb[0]*(np.pi/180.), or_axis, gb[1], target_dir='./example_dir', gb_type='tilt')

