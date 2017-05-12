import re
import os
import sys
import glob
import argparse
import subprocess

class VaspWriter(object):
  """
  Pass object a str representation of incar file.
  method mod_car is overloaded to accept a dictionary
  of variables that you may wish to modify {npar = 8}.
  """
  def __init__(self, incar_str):
    self.incar_str =incar_str

  def mod_incar(self, **kwargs):
    for k, v in kwargs.items():
      if k=='npar':
        self.incar_str = re.sub('NPAR\s+=\s+[0-9]+', 'NPAR = {npar}'.format(npar=v), self.incar_str)
      if k=='kpar':
        self.incar_str =  re.sub('KPAR\s+=\s+[0-9]+', 'KPAR = {kpar}'.format(kpar=v), self.incar_str)
      if k=='ediffg':
        self.ediff_g = 'EDIFFG = {ediffg}'.format(ediffg=v)

  def write_incar(self):
    """
    """
    with open('INCAR','w') as f:
      print >> f, self.incar_str

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--submit', action='store_true')
args = parser.parse_args()

jobs = glob.glob('b1*')
jobs = filter(lambda x: os.path.isdir(x), jobs)

scratch = os.getcwd()
for job in jobs[1:]:
  os.chdir(job)
  print 'python', os.getcwd()
  job = subprocess.Popen('cp /home/mmm0007/pymodules/templates/vasp.sh ./'.split())
  job.wait()
  with open('INCAR','r') as f:
    incar = f.read()
  with open('INCAR','w') as f:
    print >> f, re.sub('NPAR\s+=\s+4', 'NPAR   = 24', incar)
  if args.submit:
    job = subprocess.Popen('qsub vasp.sh'.split())
    job.wait()
  os.chdir(scratch)
