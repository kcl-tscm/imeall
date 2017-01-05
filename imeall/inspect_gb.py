import os
import sys
import json
import glob
import argparse
from   peewee    import *
from   quippy    import Atoms
from   datetime  import datetime, timedelta
from   models   import GBAnalysis, PotentialParameters
from   gb_models import GrainBoundary, SubGrainBoundary


GRAIN_DATABASE = "/home/lambert/pymodules/imeall/imeall/grain_boundaries/"
DATABASE       = "/home/lambert/pymodules/imeall/imeall/gb_database.db"
database       = SqliteDatabase(DATABASE)


class Inspector(object):
  """
  This class inspects the SQLdatabase to list GB energies according to potential,
  orientation axis etc., it highlights unconverged calculations. 
  """
  def __init__(self):
    pass

  def list_gb(self, potential="PotBH.xml", or_axis="001", print_unconverged=True):
    """
    list energies, and convergence for gbs for potential and orientation axis. 
    """
    pot_param     = PotentialParameters()
    ener_per_atom = pot_param.gs_ener_per_atom()
    #serialize orientation vector:
    oraxis        = ','.join([c for c in or_axis])
    gbs           = GrainBoundary.select().where(GrainBoundary.orientation_axis==oraxis).order_by(GrainBoundary.angle)
    unconverged   = []
    for gb in gbs:
      print "OR axis: {or_axis}, Angle: {angle}".format(or_axis=gb.orientation_axis, angle=round(gb.angle*(180.0/3.14159),2))
      subgbs = (gb.subgrains.select(GrainBoundary, SubGrainBoundary)
                  .where(SubGrainBoundary.potential==potential)
                  .join(GrainBoundary)
                  .order_by(SubGrainBoundary.E_gb)
                  .dicts())
      if len(subgbs) > 0:
        subgbs = [(16.02*(subgb['E_gb']-float(subgb['n_at']*ener_per_atom['PotBH.xml']))/(2.0*subgb['area']), subgb) for subgb in subgbs]
        subgbs.sort(key = lambda x: x[0])
        for subgb in subgbs:
          print "\t", subgb[1]['gbid'], subgb[0], subgb[1]['converged']
          if not subgb[1]['converged']:
            unconverged.append((subgb[1]['gbid'], subgb[1]['path']))
      else:
        print "No Subgrains."

    print "{} unconverged subgrains".format(len(unconverged))
    if print_unconverged:
      for unconv in unconverged:
        print unconv[1]

if __name__=="__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("-p", "--potential", help="potential to pull from database.", default="PotBH.xml")
  parser.add_argument("-o", "--or_axis",   help="orientation axis to pull from database.", default="001")
  parser.add_argument("-c", "--converged", help="print list of unconverged grains.", action="store_true")
  args   = parser.parse_args()
  inspector = Inspector()
  inspector.list_gb(potential=args.potential, or_axis=args.or_axis, print_unconverged=args.converged)
