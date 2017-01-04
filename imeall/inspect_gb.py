import os
import sys
import json
import glob
import argparse
from   peewee    import *
from   quippy    import Atoms
from   models    import GBAnalysis
from   datetime  import datetime, timedelta
from   gb_models import GrainBoundary, SubGrainBoundary


GRAIN_DATABASE = "/home/lambert/pymodules/imeall/imeall/grain_boundaries/"
DATABASE       = "/home/lambert/pymodules/imeall/imeall/gb_database.db"
database       = SqliteDatabase(DATABASE)

#parser  = argparse.ArgumentParser()
#parser.add_argument("-o", "--orientation_axis", default="110")
#args = parser.parse_args()

min_ens = (GrainBoundary
             .select(GrainBoundary, SubGrainBoundary)
             .join(SubGrainBoundary)
             .where(SubGrainBoundary.potential=='PotBH.xml')
             .group_by(SubGrainBoundary.canonical_grain)
             .order_by(GrainBoundary.angle)
             #.having(SubGrainBoundary.E_gb == fn.Min(SubGrainBoundary.E_gb))
             .dicts())

for minen in min_ens:
  print minen['gb_type'], map(int, minen['orientation_axis'].split(',')), minen['angle'], minen['E_gb']

