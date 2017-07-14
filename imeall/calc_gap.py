import os
import sys
import json
import glob
import logging
import argparse
import numpy as np
  
from peewee   import *
from datetime import datetime, timedelta
from models       import GBAnalysis, PotentialParameters
from gb_models    import GrainBoundary, SubGrainBoundary
from quippy       import Atoms, set_fortran_indexing, Potential, AtomsReader
from ase.optimize import FIRE
from ase.constraints import UnitCellFilter, StrainFilter

oraxis = '0,0,1'
pot_param     = PotentialParameters()
ener_per_atom = pot_param.gs_ener_per_atom()
selected_grains = GrainBoundary.select().where(GrainBoundary.orientation_axis==oraxis).where(GrainBoundary.boundary_plane != oraxis)

f = open('./locenviron/gap_energies.dat','a')

for gb in selected_grains.order_by(GrainBoundary.angle)[2:]:
  subgbs = (gb.subgrains.select(GrainBoundary, SubGrainBoundary)
              .where(SubGrainBoundary.potential=='PotBH.xml')
              .join(GrainBoundary).dicts())
  subgbs = [(16.02*(subgb['E_gb']-float(subgb['n_at']*ener_per_atom['PotBH.xml']))/(2.0*subgb['area']), subgb) for subgb in subgbs]
  subgbs.sort(key = lambda x: x[0])
  try:
    print subgbs[0][1]['path']
    continue
    
    target_dir = os.path.join('./grain_boundaries', subgbs[0][1]['path'])
    struct_file = os.path.join(target_dir, subgbs[0][1]['gbid']) + '_traj.xyz'
    print struct_file
    ats = AtomsReader(struct_file)[-1]
    pot = Potential('IP GAP', param_filename='gp33b.xml')
    ats.set_calculator(pot)
    print subgbs[0][1]['n_at'], subgbs[0][1]['area']
    strain_mask = [0,0,1,0,0,0]
    ucf = UnitCellFilter(ats, strain_mask)
    opt = FIRE(ucf)
    FORCE_TOL = 0.1
    opt.run(fmax=FORCE_TOL)
    gap_en = ats.get_potential_energy()
    print gap_en
    print round(gb.angle*(180.0/3.14159),3), round(subgbs[0][0],3), 16.02*(gap_en-float(subgbs[0][1]['n_at']*ener_per_atom['gp33b.xml']))/(2.0*subgbs[0][1]['area'])
    print >> f, round(gb.angle*(180.0/3.14159),3), round(subgbs[0][0],3), 16.02*(gap_en-float(subgbs[0][1]['n_at']*ener_per_atom['gp33b.xml']))/(2.0*subgbs[0][1]['area'])
    ats.write('./locenviron/{}.xyz'.format(subgbs[0][1]['gbid']))
  except IndexError:
    print '\t', round(gb.angle*(180.0/3.14159),3), subgbs


