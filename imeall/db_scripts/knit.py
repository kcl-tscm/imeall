import os
import sys
import glob
from quippy import AtomsReader, AtomsWriter


#Read in a collection of traj output files
#and dump them into a single traj file.

traj_files = sorted(glob.glob('[0-9]*.traj.xyz'))
knit_traj = AtomsWriter('full_traj.xyz') 

for traj_file in traj_files:
  ats = AtomsReader(traj_file)
  print len(ats)
  for at in ats:
    mom = [3.0 for x in range(len(at))]
    at.set_initial_magnetic_moments(mom)
    at.add_property('mhm1', at.properties['modified_hybrid_mark_1'])
    for i, mark in enumerate(at.properties['modified_hybrid_mark_1']):
      if mark> 0:
        print i, mark
    knit_traj.write(at)

  
  
  
