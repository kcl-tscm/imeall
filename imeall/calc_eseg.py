import sys
from ase import Atoms as aseAtoms
from ase.io.vasp import write_vasp
from ase.calculators.vasp import Vasp
from pyspglib import spglib
from quippy import Atoms, farray, frange, set_fortran_indexing
from calc_inter_dist  import decorate_interface

import json
import argparse
import numpy as np

sys.settrace

def nearest_to_unique(at, unique_sites):
  """Given an :class:`ase.Atom` object and an array of site vectors finds
  whether there is a correspondence between the atom position
  and one of the sites.

  Args:
    at(:class:`ase.Atom`): atom object. 
    unique_sites(list): list of position vectors.

  Returns:
    bool: Whether or not atom position is identical to a site
    in the unique_sites list.
  """
  equiv_site = False
  for site in unique_sites:
    if np.allclose(at.position, np.array(site[:3])):
      equiv_site = True
  return equiv_site


def gen_interface():
  """Selects an interfacial region of a bicrystal 
  based on common neighbour analysis. The width of the interfacial region
  is equal to 2*(gb_max-gb_min) where gb_max is the z-coordinate of the highest
  non-bcc atom, and gb_min is the lowest non-bcc atom.

  The method creates a file `interface.xyz` in the working directory,
  with the interface centered in a unit cell with 1 angstrom vacuum 
  on each side. 

  Returns:
    :class:`ase.Atoms`: Atoms object of the interfacial slab in same 
    coordinates as original bicrystal.
  """
  #output.xyz must have structure_type property attached.
  ats = Atoms('output.xyz')
  cell_midpoint = ats.get_cell()[2,2]/2.0
  #select non-BCC sites are 0 otherwise 3.
  struct_type = np.array(ats.properties['structure_type'])
  struct_mask = [not struct for struct in struct_type]
  interface = ats.select(struct_mask)
  #select upper interface to decorate.
  interface = interface.select([at.position[2] > cell_midpoint for at in interface])
  z_vals = [at.position[2] for at in interface]
  z_min = min(z_vals)
  z_max = max(z_vals)
  #take slice of interface max uncoordinated with min uncoordinated.
  z_width = (z_max-z_min)/2.0
  z_center = z_width + z_min
  gb_max = z_max + 1.0*z_width
  gb_min = z_min - 1.0*z_width
  zint = ats.select([(gb_min <= at.position[2] <= gb_max) for at in ats])
  #make a copy to return
  int_ats = zint.copy()
  zint.center(vacuum=1.0, axis=2)
  zint.write('interface.xyz')
  #Write POSCAR to use interstitial site generator:
  ats = Atoms('interface.xyz')
  vasp_args=dict(xc='PBE', amix=0.01, amin=0.001, bmix=0.001, amix_mag=0.01, bmix_mag=0.001,
               kpts=[3, 3, 3], kpar=9, lreal='auto', ibrion=-1, nsw=0, nelmdl=-15, ispin=2,
               nelm=100, algo='VeryFast', npar=24, lplane=False, lwave=False, lcharg=False, istart=0,
               voskown=0, ismear=1, sigma=0.1, isym=2)
  vasp = Vasp(**vasp_args)
  vasp.initialize(ats)
  write_vasp('POSCAR', vasp.atoms_sorted, symbol_count=vasp.symbol_count, vasp5=True)
  return int_ats

if __name__=='__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('--gen_interface','-g',help='create a slab of the interfacial region.', action='store_true')
  args = parser.parse_args()

  if args.gen_interface:
    gen_interface()

  if args.decorate_interface:
    decorate_interface()

