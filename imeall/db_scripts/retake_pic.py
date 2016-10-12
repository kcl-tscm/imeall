import os
import sys
import glob
import argparse
import imeall.slabmaker.slabmaker as slabmaker

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--pattern", help="Job pattern.")
args = parser.parse_args()

scratch = os.getcwd()
pattern = args.pattern
jobdirs = ['./']

for job in jobdirs:
  os.chdir(job)
  struct_file = glob.glob(pattern+'.xyz')
  print job
  var = raw_input('Retake photo?')
  if var =='y':
    print 'retaking photo'
    slabmaker.take_pic(pattern, translate=True)
  elif var =='n':
    pass

