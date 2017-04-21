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
parser.add_argument("-gbt","--gb_type", help="Specify type of boundary twist or tilt.", required=True)
parser.add_argument("-ct","--calc_type", help="Potential used for calculation.", required=True) 
parser.add_argument("-p","--pattern", help="Job pattern to select grainboundaries.", default="001")
parser.add_argument("-d","--delay", help="Time delay between job submissions.", type=int, default=60)
parser.add_argument("-t","--tranchsize", help="Tranchsize.", type=int, default=2000)
parser.add_argument("-q","--queue", help="Name of queue on machine.", default="LowMemShortterm.q") 
parser.add_argument("-j","--jobfile", help="Job file if this is specified the jobs will only be run\
                                            if they are present in the file specified. Each gbid \
                                            should be on a different line.", default="") 
args = parser.parse_args()

#completed_jobs = ['00195301120','00181701140','0017575790','00171501160','0016193350','00149556130','00111421100']
#[jobdirs.remove(x) for x in completed_jobs]

def chunker(seq, size):
  return (seq[pos:pos + size] for pos in xrange(0, len(seq), size))

if args.gb_type=="tilt":
  if args.jobfile=="":
    jobdirs   = glob.glob('{dir_pattern}*'.format(dir_pattern=args.pattern))
    jobdirs   = filter(os.path.isdir, jobdirs)
  else:
    with open(args.jobfile, 'r') as f:
      jobdirs = f.read().split()
  print jobdirs
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
  jobdirs   = filter(os.path.isdir, jobdirs[1:])
  jdirs     = []
  job_index = 0
  for job in jobdirs:
    for rc in np.arange(1.1, 2.3, 0.1):
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
for job_tract in chunker(jobdirs, args.tranchsize):
  for job in job_tract:
    os.chdir(os.path.join(scratch, job[0]))
    print 'Current Working Directory', os.getcwd()
    if parallel:
      print 'Parallel Job'
      pbs_str = open('/users/k1511981/pymodules/templates/calc_rundyn.pbs', 'r').read()
      gb_args = '-ct {calc_type} -rc {rc} -i_v {i_v} -i_bxv {i_bxv} -gbt {gb_type}'.format(rc=job[1], i_v=job[2], i_bxv=job[3], 
      calc_type = args.calc_type, gb_type=args.gb_type)
      pbs_str = pbs_str.format(jname='fe'+job[0][:8], time=jtime, queue=args.queue, gb_args=gb_args)
      with open('gb.pbs', 'w') as pbs_file:
        print >> pbs_file, pbs_str
      qsub_args = 'qsub gb.pbs'
      print >> log, job, gb_args, job[4]
      job     = subprocess.Popen(qsub_args.split())
      job.wait()
    else:
      print 'Running in Serial'
      gb_args = '-ct {calc_type} -rc {rc} -i_v {i_v} -i_bxv {i_bxv} -gbt {gb_type}'.format(rc=job[1], i_v=job[2], i_bxv=job[3],
                gb_type=args.gb_type, calc_type=args.calc_type)
      gb_args = "python /Users/lambert/pymodules/imeall/imeall/run_dyn.py {}".format(gb_args)
      print job, gb_args
      job = subprocess.Popen(gb_args.split())
      job.wait()
  time.sleep(args.delay)
  os.chdir(scratch)
  log.flush()
log.close()
