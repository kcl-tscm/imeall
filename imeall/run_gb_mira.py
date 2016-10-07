#!/usr/bin/env python
import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)

import os
import sys
import glob
import socket
import subprocess
import numpy as np

from bgqtools import (get_bootable_blocks, boot_blocks, block_corner_iter,
                      set_unbuffered_stdout)

set_unbuffered_stdout()

acct    = 'SiO2_Fracture'
time    = 60
queue   = 'default'
mapping = 'ABCDET'
scratch = os.getcwd()
vasp    = '/projects/SiO2_Fracture/iron/vasp.bgq'
envargs = '--envs RUNJOB_MAPPING=%s --envs MPIRUN_ENABLE_TTY_REPORTING=0' % mapping
npj     = 1 # nodes per job
ppn     = 8 # MPI tasks per node

hostname = socket.gethostname()
print 'Hostname: %s' % hostname

jobdirs   = glob.glob('./grain_boundaries/110*')
jobdirs = filter(os.path.isdir, jobdirs)
print 'jobdirs = %s' % jobdirs

for job in jobdirs:
  for rc in np.arange(1.7, 1.79, 0.1):
    for i in np.arange(0.0, 0.51, 0.1):
      for j in np.arange(0.0, 0.51, 0.1):
        jobdirs.append((job, rc,i,j))

njobs = len(jobdirs)
nodes = npj*njobs

if 'COBALT_PARTSIZE' not in os.environ:
    print 'Not running under control of cobalt. Launching qsub...'
    qsub_args = 'qsub -A %s -n %d -t %d -q %s --mode script --disable_preboot %s' % (acct, nodes, time, queue, ' '.join(sys.argv))
    print qsub_args
    os.system(qsub_args)
    sys.exit(1)

partsize  = int(os.environ['COBALT_PARTSIZE'])
partition = os.environ['COBALT_PARTNAME']
jobid     = int(os.environ['COBALT_JOBID'])

print 'Nodes per job: %d' % npj
print 'MPI tasks per node: %d' % ppn
print 'Number of sub-block jobs: %d' % njobs
print 'Number of nodes: %d' % nodes

blocks = get_bootable_blocks(partition, nodes)
print 'Available blocks: %s' % blocks
boot_blocks(blocks)
# start sub-block jobs with background runjob helper processes
jobs = []
logs = []

# Lets make each job a tuple struct like
# (jobdir, bxv, v, rc).
for job, (block, corner, shape) in zip(jobdirs, block_corner_iter(blocks, npj)):
  print job, (block, corner, shape)
  os.chdir(os.path.join(scratch, job[0]))
  log = open('gbrelax.out', 'w')
  locargs = '--block %s --corner %s --shape %s' % (block, corner, shape)
# runjob_args = ('python %s -n %d -p %d %s' % (locargs, npj*ppn, ppn, envargs)).split()
  pyargs  = 'python -rc {rc} -i_v {i_v} -i_bxv {i-bxv} '.format(rc=job[1], i_v=job[2], i_bxv=job[3])
  runjob_args = ('runjob %s -n %d -p %d %s : %s'%(locargs, 1, 1, envargs, pyargs)).split()
  print ' '.join(runjob_args)
  jobs.append(subprocess.Popen(runjob_args, stdout=log))
  logs.append(log)
    
# wait for all background jobs to finish, then flush their logs
for (job, log) in zip(jobs, logs):
  job.wait()
  log.flush()
