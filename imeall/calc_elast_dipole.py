import os
import subprocess

hydr_py = '/Users/lambert/pymodules/imeall/imeall/hydrogenate_gb.py'
jobs = [('bulk',1.00),('bulk_compress', 0.98),('bulk_strained', 1.02), ('gb',1.00),('gb_compress', 0.98),('gb_strained', 1.02)]
for job in jobs:
  if 'gb' in job[0]:
    launch = subprocess.Popen('python {hydr_py} -s {strain} -d {dirname}'.format(hydr_py=hydr_py, strain=job[1], dirname=job[0]).split())
  else:
    launch = subprocess.Popen('python {hydr_py} -s {strain} -d {dirname} -b'.format(hydr_py=hydr_py, strain=job[1], dirname=job[0]).split())
  launch.wait()




