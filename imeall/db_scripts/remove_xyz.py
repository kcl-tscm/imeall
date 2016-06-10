import os
import sys
from models import GBMaintenance 

gbm = GBMaintenance()
for thing in os.listdir('./'):
  if os.path.isdir(thing) and thing[:len(sys.argv[1])] == sys.argv[1]:
    gbm.remove_xyz(thing, dryrun=False)
    
