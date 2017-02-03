import os
import bz2
import glob
import shutil
import subprocess
import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-p","--prefix",help="prefix", required=True)
args = parser.parse_args()

def _comp_sub(filename):
  print 'Compressing {0}'.format(filename)
  job_args = "gzip -6 {0}".format(filename)
  job = subprocess.Popen(job_args.split())
  job.wait()

jobdirs = filter(os.path.isdir, glob.glob("{0}*".format(args.prefix)))
scratch = os.getcwd()
for job in jobdirs:
  print job
  os.chdir(job)
  #traj_files = glob.glob("*.traj.xyz")
  traj_files = glob.glob("crack_traj.xyz")
  print traj_files
  for traj in traj_files:
    _comp_sub(traj)
  os.chdir(scratch)
