import imeall.slabmaker.slabmaker as slabmaker
import os
import sys

jobdirs = []
for thing in os.listdir('./'):
  if os.path.isdir(thing) and thing[:3]=='110':
    jobdirs.append(thing)

for job in jobdirs[5:]:
    print job
    var = raw_input('Retake photo?')
    if var =='y':
      slabmaker.take_pic(os.path.join(job,job), translate=True)
      print 'retaking photo'
    elif var =='n':
      pass

