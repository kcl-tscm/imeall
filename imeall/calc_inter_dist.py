import json
import spglib
import numpy as np

from ase.geometry import geometry
from pylada_defects_2 import get_unique_ints, Voronoi
from quippy import Atoms
from quippy import set_fortran_indexing

set_fortran_indexing(False)

ats = Atoms('interface.xyz')
dataset = spglib.get_symmetry_dataset(ats, symprec=1e-5)

with open('unique_lattice_sites.json','w') as f:
  json.dump([list(ats[site_num].position) for site_num in np.unique(dataset['equivalent_atoms'])], f)

unique_atoms = []
for at in ats:
  unique_atoms.append(at.position)
voronoi = Voronoi(limits=tuple(np.diag(ats.cell)), periodic=(True, True, False))
cntr = voronoi.compute_voronoi(unique_atoms)
ints_list = []
for site_num in np.unique(dataset['equivalent_atoms']):
  for vert in voronoi.get_vertices(site_num, cntr):
    ints_list.append(vert.tolist())
for unique in ints_list:
    ats.add_atoms(unique, 1)
for i in range(len(ats)):
  ats.id[i] = i
#remove voronoi duplicates
print 'Fe_H atoms', len(ats)
ats.wrap()
geometry.get_duplicate_atoms(ats, cutoff=0.2, delete=True)
print 'Fe_H atoms remove duplicates', len(ats)
#select unique hydrogens
for i in range(len(ats)):
  ats.id[i] = i

ints_list = [at.position for at in ats if at.number==1]
with open('unique_h_sites.json', 'w') as f:
  json.dump([list(u)for u in ints_list], f)
ats.write('hydrogenated_grain.xyz')

