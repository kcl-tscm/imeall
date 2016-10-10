import os
import re
import glob 
import json
import argparse
import subprocess
import ase.units as units

def create_json(dir_pattern=''):
  """
  :method:`create_json` pulls data from a vasp OUTCAR file.
  """
  jobdirs = glob.glob(dir_pattern)
  jobdirs = filter(os.path.isdir, jobdirs)
  for job in jobdirs:
    os.chdir(job)
    print job
    job = subprocess.Popen('python /Users/lambert/pymodules/imeall/imeall/db_scripts/read_vasp.py'.split())
    job.wait()
    os.chdir('../')

def plot_ener(eners, calc_pattern='', A=181.221, reverse=False):
  """
  :method:`plot_ener` plot energy.
  """
  ref_en = eners[0][1]
  with open('{0}.dat'.format(calc_pattern), 'w') as f: 
    for ener in eners:
      print >> f, ener[0], ((ener[1]-ref_en)/A)*(units.m**2/units.J)
    if reverse:
      eners.reverse()
      for ener in eners[1:]:
        print >> f, 1.-float(ener[0]), ((ener[1]-ref_en)/A)*(units.m**2/units.J)
  return

def grab_eners(dir_pattern=''):
  """
  :method:`grab_ener` plot energy.
  """
  jobdirs = glob.glob(dir_pattern)
  jobdirs = filter(os.path.isdir, jobdirs)
  eners = []
  re_disp = re.compile(r'[0-9\.]+',re.S)
  for job in jobdirs:
    print job
    #disp = re_disp.findall(job)
    os.chdir(job)
# Pull energy from the json file.
    try:
      with open('subgb.json','r') as f:
        gb_dict = json.load(f)
      eners.append((job[-3:], gb_dict['E_gb']))
    except:
      print 'Unable to read subgb', job
      pass
    os.chdir('../')
  return eners

if __name__=='__main__':
  parser  = argparse.ArgumentParser()
  parser.add_argument("-A", "--area",    help="Area of the grain boundary.", type=float, default=181.221)
  parser.add_argument("-p", "--pattern", help="glob pattern for directories containing gamma surface information.", required=True)
  parser.add_argument("-r", "--reverse", help="If energetics are symmetric can create the full plot by reversing.", action='store_true')
  args    = parser.parse_args()
  area    = args.area
  pattern = args.pattern
  reverse = args.reverse

  create_json(dir_pattern='{0}*'.format(pattern))
  eners         = grab_eners('{0}*'.format(pattern))
  calc_pattern  = '{0}'.format(pattern)
  print 'Writing data to ', calc_pattern+'.dat'
  plot_ener(eners=eners, calc_pattern=calc_pattern, A=area, reverse=reverse)
