import os
import glob
import argparse
from quippy import Atoms
from quippy.io import AtomsWriter, AtomsReader

parser = argparse.ArgumentParser()

parser.add_argument("-p", "--pattern", required=True)
parser.add_argument("-j", "--jobfile")
pattern = args.pattern

jobs = glob.glob(pattern)
out = AtomsWriter('traj.xyz')
for job  in jobs:
  struct_file = os.path.join(job,'feb{0}.xyz'.format(job[-3:]))
  ats = AtomsReader(struct_file)[-1]
  out.write(ats)
