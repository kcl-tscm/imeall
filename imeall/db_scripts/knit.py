import os
import sys
import glob
from   quippy import AtomsReader, AtomsWriter, set_fortran_indexing


#Read in a collection of traj output files
#and dump them into a single traj file.
#list of properties that can be deleted
del_prop_lit = ['hybrid_vec', 'hybrid_1', 'weight_region1', 'weight_region1_1']
traj_files = sorted(glob.glob('[0-9]*.traj.xyz'))
knit_traj = AtomsWriter('full_traj.xyz') 

for traj_file in traj_files:
  ats = AtomsReader(traj_file)
  print len(ats)
  for at in ats:
    mom = [3.0 for x in range(len(at))]
    at.set_initial_magnetic_moments(mom)
    at.add_property('mhm1', at.properties['modified_hybrid_mark_1'])
    at.add_property('edgex', 0.0, overwrite=True)
    at.add_property('edgey', 0.0, overwrite=True)
    at.add_property('screw', 0.0, overwrite=True)
    knit_traj.write(at)

  
  
  
