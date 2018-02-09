import os
import sys
from imeall.models import GBMaintenance

gbm = GBMaintenance()
for thing in os.listdir('./'):
    if os.path.isdir(thing) and thing[:len(sys.argv[1])] == sys.argv[1]:
        gbm.fix_json(thing)
