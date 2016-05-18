import os
import sys
from imeall.models import GBMaintenance 

gbm = GBMaintenance()
for thing in os.listdir('./'):
  if os.path.isdir(thing) and thing[:3] == sys.argv[1]:
    gbm.remove_eo_files(thing)
