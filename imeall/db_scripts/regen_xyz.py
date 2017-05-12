#!/usr/bin/env python
import logging
logging.basicConfig(format='%(asctime)s %(message)s', level=logging.INFO)
import os
import sys
import glob
import json
import time
import socket
import argparse
import subprocess
import numpy as np
from   datetime import datetime
from models import PotentialParameters

npj     = 1 # nodes per job

scratch = os.getcwd()

parser = argparse.ArgumentParser()
parser.add_argument("-i","--inputfile", help="Job file if this is specified the jobs will only be run\
                                            if they are present in the file specified. Each gbid \
                                            should be on a different line.", default="") 
parser.add_argument("-gbt","--gb_type", help="Specify type of boundary twist or tilt.", required=True)
parser.add_argument("-t", "--time", help="job time %H:%M:%s", default="2:00:00")
parser.add_argument("-q","--queue", help="Name of queue on machine.", default="LowMemShortterm.q") 

potparams = PotentialParameters()

args = parser.parse_args()
with open(args.inputfile, 'r') as f:
  jobdirs = f.read().split()

scratch = os.getcwd()
for job in jobdirs:
  os.chdir(os.path.join(scratch, job))
  with open('subgb.json','r') as f:
    json_dict = json.load(f)
  gbid  = json_dict['gbid']
  rc  = json_dict['rcut']
  i_v = json_dict['rbt'][0]
  i_bxv = json_dict['rbt'][1]
  paramfile = json_dict['param_file']
  canonical_gb_path = os.path.join(scratch, '/'.join(job.split('/')[0:3]))
  print canonical_gb_path
  os.chdir(canonical_gb_path)
  calc_type = potparams.potdir_dict()[paramfile]
  print calc_type
  pbs_str = open('/users/k1511981/pymodules/templates/calc_rundyn.pbs', 'r').read()
  gb_args = '-ct {calc_type} -rc {rc} -i_v {i_v} -i_bxv {i_bxv} -gbt {gb_type}'.format(rc=rc, i_v=i_v, i_bxv=i_bxv,
              calc_type = calc_type, gb_type=args.gb_type)
  pbs_str = pbs_str.format(jname='fe'+gbid[:8], time=args.time, queue=args.queue, gb_args=gb_args)
  with open('gb.pbs', 'w') as pbs_file:
    print >> pbs_file, pbs_str
  print gb_args
  qsub_args = 'qsub gb.pbs'
  job     = subprocess.Popen(qsub_args.split())
  job.wait()
  os.chdir(scratch)
