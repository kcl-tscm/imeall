import os
import sys
import re
import time
from collections import deque
import argparse

class Job(object):
  def __init__(self, job_que=deque()):
    self.pbs_file = ''
    self.job_dir  = ''
    self.job_id   = ''
    self.job_que  = job_que

  def gen_pbs(self, time='8:00:00', queue='serial.q'):
    ''' 
      Generates a pbs file from the template stored in pbs_str. Writes 
    '''
    pbs_str = open('/users/k1511981/pymodules/templates/calc_ada.pbs','r').read()
    pbs_str = pbs_str.format(jname='fe'+self.name, xyz_file='{0}.xyz'.format(self.name),
                             time=time, queue=queue)
    print os.path.join(self.subgrain_dir, 'fe{0}.pbs'.format(self.name))
    with open(os.path.join(self.subgrain_dir, 'fe{0}.pbs'.format(self.name)) ,'w') as pbs_file:
      print >> pbs_file, pbs_str

  def sub_pbs(self, job_dir, exclude, suffix='v6bxv2z', regex=None, calc_type=None):
    ''' 
    Given an explicit suffix, or a regex this routine recurses through
    the directory structure and submits any pbs files that 
    directories that (we mightn't want for instance DFT on Ada
    or EAM on Mira.)
    Useful submission patterns include:
    REGEX:
      submit all pbs files with a rigid body translation 
      and any atom deletion criterion hydrogen concentration etc.: 
        tv[.0-9]+bxv[.0-9]+_.*?
      submit all sub-pbs files with a rigid body translation
      and different atom _deletion criterion: 
        tv[.0-9]+bxv[.0-9]+_d[.0-9]+
      all translations with a specific deletion criterion
      in this case 2.3 A:
        tv[.0-9]+bxv[.0-9]+_d2.3
      etc.
    SUFFIX:
      submit all super cells: 
        v6bxv2z
      '''
    lst = os.listdir(job_dir)
    for dir in lst:
      dir = os.path.join(job_dir, dir)
      if regex == None:
        if os.path.isdir(dir) and dir != 'DFT':
          print dir
          self.sub_pbs(dir, exclude, suffix=suffix, regex=regex)
        elif dir.split('_')[-1] == suffix:
          pbs_dir = os.path.join(sub_dir, dir)
        else:
          pass
      else:
        if os.path.isdir(dir) and dir.split('/')[-1] not in exclude: 
          self.sub_pbs(dir,exclude, suffix=suffix, regex=regex)
        elif regex.match(dir):
          try:
            dir  = '/'.join(dir.split('/')[:-1])
            name = dir.split('/')[-1]
            self.job_que.append([dir, name])
          except:
            print 'Job Submit Failed'
        else:
          pass

if __name__=='__main__':
#String for orientation axis e.g. '110'
  jobdirs = []
  que = deque()
  job = Job(job_que=que)
  parser = argparse.ArgumentParser()
  parser.add_argument("-p", "--prefix", help="Choose directories to act on", required=True)
  parser.add_argument("-ct", "--calc_type", help="Choose Calculation Type to submit", required=True)
  parser.add_argument("-rc", "--rcut", help="Choose cutoff radius to submit", required=True)
  parser.add_argument("-ht", "--hold_time", help="Set wait time between tranches of submitted jobs (seconds)", default=60, type=int)

  args = parser.parse_args()
  prefix     = args.prefix
  calc_type  = args.calc_type
  hold_time  = args.hold_time
  rcut       = args.rcut

  for thing in os.listdir('./'):
    if os.path.isdir(thing) and thing[:len(prefix)]==prefix:
      jobdirs.append(thing)

  exclude = ['EAM_Mish', 'DFT', 'EAM_Men','EAM_Ack','EAM_Dud', 'EAM']
  exclude.remove(calc_type)

  #suffix_reg = re.compile(r'.*?v6bxv2_tv0.0+bxv0.0+_d2.3z.pbs')
  #suffix_reg = re.compile(r'.*?v6bxv2_tv0.0+bxv0.0+_d{0}z.pbs'.format(rcut))
  suffix_reg = re.compile(r'.*?v6bxv2_tv[0-9.]+bxv[0-9.]+_d{0}z.pbs'.format(rcut))
  for thing in jobdirs[:]:
    job.sub_pbs(thing, exclude, regex=suffix_reg, calc_type=calc_type)

  counter = 0
  while que:
    for i in range(32):
      try:
        que_empty = False
        a = que.popleft()
        os.system("cd {0}; qsub fe{1}.pbs".format(a[0], a[1]))
        counter += 1
      except IndexError:
        que_empty = True
        print 'job queu empty!'
      if que_empty: 
        break
    if que: 
      print 'Submitted', counter, 'jobs'
      print 'taking a break'
      time.sleep(hold_time)
