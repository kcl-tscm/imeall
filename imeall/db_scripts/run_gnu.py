import os
import sys

jobdirs = []
for thing in os.listdir('./'):
  if os.path.isdir(thing) and thing[:3]=='111':
    jobdirs.append(thing)
for job_dir in jobdirs:
  os.system('cd {0}; gnuplot {1}'.format(job_dir, 'plot.gnu'))

