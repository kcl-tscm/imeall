import os
import subprocess

example_dirs = ['boundary_interstitials', 'checkout_structures', 'dissolution_test', 
'elastic_dipole_test','gen_grainboundaries', 'query_test']

scratch = os.getcwd()

for dir_ in example_dirs:
  os.chdir(dir_)
  print dir_
  job = subprocess.Popen('sh run.sh'.split())
  job.wait()
  os.chdir(scratch)
