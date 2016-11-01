import os
import subprocess

jobs = [('bulk',1.00),('bulk_compress', 0.98),('bulk_strained', 1.02), ('gb',1.00),('gb_compress', 0.98),('gb_strained', 1.02)]
for job in jobs:
  if 'gb' in job[0]:
    launch = subprocess.Popen('python hydrogenate_gb.py -s {strain} -d {dirname}'.format(strain=job[1], dirname=job[0]).split())
  else:
    launch = subprocess.Popen('python hydrogenate_gb.py -s {strain} -d {dirname} -b'.format(strain=job[1], dirname=job[0]).split())
  launch.wait()




