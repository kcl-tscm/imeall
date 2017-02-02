#!/usr/bin/env python
import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
import os
import sys
import glob
import time
import socket
import argparse
import subprocess
import numpy as np
from   datetime import datetime

npj     = 1 # nodes per job
scratch = os.getcwd()

parser = argparse.ArgumentParser()
parser.add_argument("-gbt", "--gb_type", help="Specify type of boundary twist or tilt.", default="tilt")
parser.add_argument("-p",   "--pattern", help="Job pattern to select grainboundaries.", default="001")
parser.add_argument("-d",   "--delay",   help="Time delay between job submissions.", type=int, default=60)
args = parser.parse_args()

def chunker(seq, size):
  return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))

if args.gb_type=="tilt":
  jobdirs   = glob.glob('{}*'.format(args.pattern))
  jobdirs   = filter(os.path.isdir, jobdirs)
  jdirs     = []
  job_index = 0
  for job in jobdirs:
    for rc in np.arange(1.4, 2.1, 0.1):
      for i in np.arange(0.0, 0.50, 0.1):
        for j in np.arange(0.0, 0.50, 0.1):
          job_index += 1
          jdirs.append((job, rc, i, j, job_index))
  jobdirs = filter(lambda x: os.path.isdir(x[0]), jdirs)
elif args.gb_type=="twist":
  jobdirs   = glob.glob('{}*'.format(args.pattern))
  jobdirs   = filter(os.path.isdir, jobdirs)
  print jobdirs
  jdirs     = []
  job_index = 0
  for job in jobdirs:
    for rc in np.arange(1.4, 2.3, 0.1):
      job_index += 1
      jdirs.append((job, rc, 0.0, 0.0, job_index))
  jobdirs = filter(lambda x: os.path.isdir(x[0]), jdirs)
else:
  print "Only tilt or twist boundaries valid."
  sys.exit()

njobs  = len(jobdirs)
jtime  = '2:00:00'
#test for parallel environ qsub
which_qsub = subprocess.call("which qsub".split())
if which_qsub == 0:
  parallel = True
else: 
  parallel = False

log  = open('gbcalclog.out', 'a')
for job_tract in chunker(jobdirs, 72):
  for job in job_tract:
    os.chdir(os.path.join(scratch, job[0]))
    print 'Current Working Directory', os.getcwd()
    if parallel:
      pbs_str = open('/users/k1511981/pymodules/templates/calc_rundyn.pbs', 'r').read()
      gb_args = '-rc {rc} -i_v {i_v} -i_bxv {i_bxv} -gbt {gb_type}'.format(rc=job[1], i_v=job[2], i_bxv=job[3], gb_type=args.gb_type)
      pbs_str = pbs_str.format(jname='fe'+job[0][:8], time=jtime, queue='LowMemShortterm.q', gb_args=gb_args)
      with open('gb.pbs', 'w') as pbs_file:
        print >> pbs_file, pbs_str
      qsub_args = 'qsub gb.pbs'
      print >> log, job, gb_args, job[4]
      job     = subprocess.Popen(qsub_args.split())
      job.wait()
    else:
      gb_args = '-rc {rc} -i_v {i_v} -i_bxv {i_bxv} -gbt {gb_type}'.format(rc=job[1], i_v=job[2], i_bxv=job[3], gb_type=args.gb_type)
      gb_args = "python /Users/lambert/pymodules/imeall/imeall/run_dyn.py {}".format(gb_args)
      print job, gb_args
      job = subprocess.Popen(gb_args.split())
      job.wait()
  time.sleep(args.delay)
  os.chdir(scratch)
  now       = datetime.now()
  log.flush()
log.close()
