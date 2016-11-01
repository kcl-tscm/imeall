#!/usr/bin/env python
import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

import os
import sys
import glob
import socket
import subprocess
import numpy as np
import time

from datetime import datetime

npj     = 1 # nodes per job
scratch = os.getcwd()

def chunker(seq, size):
  return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))

jobdirs = glob.glob('001*')
#print 'jobdirs = %s' % jobdirs
jobdirs = filter(os.path.isdir, jobdirs)
jdirs = []
job_index = 0
for job in jobdirs:
  for rc in np.arange(1.6, 2.31, 0.1):
    for i in np.arange(0.0, 0.51, 0.1):
      for j in np.arange(0.0, 0.51, 0.1):
        job_index += 1
        jdirs.append((job, rc, i, j, job_index))

jobdirs = filter(lambda x: os.path.isdir(x[0]), jdirs)
jobdirs = jobdirs[1850:]
njobs   = len(jobdirs)
nodes   = npj*njobs
jtime    = '1:00:00'
print len(jobdirs)

jobs = []
logs = []

sync_log = open('sync_log.txt','a')

for job_tract in chunker(jobdirs, 36):
  for job in job_tract:
    os.chdir(os.path.join(scratch, job[0]))
    log = open('gbcalclog.out', 'a')
    pbs_str = open('/users/k1511981/pymodules/templates/calc_rundyn.pbs', 'r').read()
    gb_args = '-rc {rc} -i_v {i_v} -i_bxv {i_bxv} '.format(rc=job[1], i_v=job[2], i_bxv=job[3])
    pbs_str = pbs_str.format(jname='fe'+job[0][:8], time=jtime, queue='smp.q', gb_args=gb_args)
    with open('gb.pbs', 'w') as pbs_file:
      print >> pbs_file, pbs_str
    qsub_args = 'qsub gb.pbs'
    print >> log, job, gb_args, job[4]
    job = subprocess.Popen(qsub_args.split())
    job.wait()
  time.sleep(1500)
#sync to local machine:
  os.chdir(scratch)
  now = datetime.now()
  print >> sync_log, now.strftime("%Y-%m-%d %H:%M:%S")
  sync_job = subprocess.Popen("rsync -auv --exclude-from rsync_exclude.txt ./ lambert@137.73.4.204:/Users/lambert/pymodules/imeall/imeall/grain_boundaries/alphaFe/001".split(), stdout=sync_log)
  sync_job.wait()
#clean up on ada so we don't go over quota.
  script_root = "/users/k1511981/pymodules/imeall/imeall/db_scripts/remove_xyz.py"
  rm_xyz      = subprocess.Popen("python {0} -p 001 -rt all ".format(script_root).split())
  rm_xyz.wait()
