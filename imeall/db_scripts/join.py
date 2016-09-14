import os
import glob
import argparse
from quippy import Atoms
from quippy.io import AtomsWriter, AtomsReader

pattern = '[0-9]*.traj.xyz'
j_files = sorted(glob.glob(pattern))
out = AtomsWriter('traj.xyz')
for j_file  in j_files:
  for at in AtomsReader(j_file):
    out.write(at)
