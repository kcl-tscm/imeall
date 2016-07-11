import os
import sys

jobdirs = []

prefix = sys.argv[1]

for thing in os.listdir('./'):
  if os.path.isdir(thing) and thing[:len(prefix)]==prefix:
    jobdirs.append(thing)
for job_dir in jobdirs:
  os.system('cd {0}; gnuplot {1}'.format(job_dir, 'plot.gnu'))

