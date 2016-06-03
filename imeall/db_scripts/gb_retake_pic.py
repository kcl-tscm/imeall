import os
import sys
import imeall.slabmaker.slabmaker as slabmaker
from   imeall.models import GBMaintenance

jobdirs = []
gb_main = GBMaintenance()
for thing in os.listdir('./'):
  if os.path.isdir(thing) and thing[:3]==sys.argv[1]:
    jobdirs.append(thing)

for job in jobdirs[8:]:
  print job
  gb_main.retake_pic(job, translate=True, toggle= True)
  
