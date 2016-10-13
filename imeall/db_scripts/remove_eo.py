import os
import sys
from   models import GBMaintenance 

gbm = GBMaintenance()
jobdirs = filter(os.path.isdir, os.listdir('./'))

for thing in jobdirs:
  gbm.remove_eo_files(thing)

