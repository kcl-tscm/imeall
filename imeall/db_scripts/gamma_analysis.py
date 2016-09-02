import os
import re
import glob 
import json
import subprocess
import ase.units as units

def create_json(dir_pattern=''):
  jobdirs = glob.glob(dir_pattern)
  for job in jobdirs:
    os.chdir(job)
    job = subprocess.Popen('python /Users/lambert/pymodules/imeall/imeall/db_scripts/read_vasp.py'.split())
    job.wait()
    os.chdir('../')

def plot_ener(eners, calc_pattern=''):
  A = 11.32629500106
  ref_en = eners[0][1]
  with open('{0}.dat'.format(calc_pattern), 'w') as f: 
    for ener in eners:
      print >> f, ener[0], ((ener[1]-ref_en)/A)*(units.m**2/units.J)
    eners.reverse()
    for ener in eners[1:]:
      print >> f, 1.-float(ener[0]), ((ener[1]-ref_en)/A)*(units.m**2/units.J)
  return

def grab_eners(dir_pattern=''):
  jobdirs = glob.glob(dir_pattern)
  eners = []
  re_disp = re.compile(r'[0-9\.]+',re.S)
  for job in jobdirs:
    disp = re_disp.findall(job)
    os.chdir(job)
    with open('subgb.json','r') as f:
      gb_dict = json.load(f)
    eners.append((disp[1], gb_dict['E_gb']))
    os.chdir('../')
  return eners

if __name__=='__main__':
  create_json(dir_pattern='b111*')
  dir_pattern = 'b111*'
  eners       = grab_eners(dir_pattern)
  calc_pattern = 'gam_110111_dft'
  print 'Writing data to', calc_pattern+'.dat'
  plot_ener(eners=eners, calc_pattern=calc_pattern)


