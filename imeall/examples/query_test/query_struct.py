#example analysis script prints angles and lowest energy for a user chosen orientation axis
#and a list of potentials, this duplicates the functionality of the underlying

import argparse
import numpy as np
from collections import OrderedDict
from imeall.models import PotentialParameters
from imeall.gb_models import database, GrainBoundary, SubGrainBoundary

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--pots', nargs='+', help='Potentials to pull energies for.', default=['Fe_Dudarev.xml', 'PotBH.xml', 'iron_mish.xml'])
parser.add_argument('-or', '--or_axis', help='Potentials to pull energies for.',default="0,0,1")
args = parser.parse_args()

database.connect()
pot_param     = PotentialParameters()
ener_per_atom = pot_param.gs_ener_per_atom()
gbs = GrainBoundary.select().where(GrainBoundary.orientation_axis==args.or_axis).where(GrainBoundary.boundary_plane != args.or_axis)
for gb in gbs.order_by(GrainBoundary.angle):
  pot_dict = OrderedDict({})
  for potential in args.pots:
    subgbs = (gb.subgrains.select(GrainBoundary, SubGrainBoundary)
                  .where(SubGrainBoundary.potential==potential)
                  .join(GrainBoundary)
                  .order_by(SubGrainBoundary.E_gb)
                  .dicts())
    #for subgb in subgbs:
    # print subgb['path']
    subgbs = [(16.02*(subgb['E_gb']-float(subgb['n_at']*ener_per_atom[potential]))/(2.0*subgb['area']), subgb) for subgb in subgbs]
    subgbs.sort(key = lambda x: x[0])
    pot_dict[potential] = subgbs[0][0]
  print '{:.3f}'.format(180.0/np.pi*gb.angle), ' '.join(['{:.3f}'.format(x) for x in pot_dict.values()])
