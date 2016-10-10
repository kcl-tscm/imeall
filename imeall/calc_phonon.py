import os
import json
import argparse
import numpy as np

from ase.io   import vasp
from ase      import Atoms as aseAtoms
from quippy   import Atoms, Potential, Atoms
from phonopy  import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from ase.lattice.cubic       import BodyCenteredCubic

def calc_band_paths(band_paths, npoints=100): 
  """
  Given a list of special points calculates the
  q-path between them.
  """
  bands = []
  nd = len(band_paths)
  print nd
  for i in range(nd - 1):
    print band_paths[i+1], band_paths[i]
    diff = (band_paths[i + 1] - band_paths[i]) / float(npoints)
    band = [band_paths[i].copy()]
    q = np.zeros(3)
    for j in range(npoints):
      q += diff
      band.append(band_paths[i] + q)
    bands.append(band)
  return bands

#CHOOSE POTENTIAL
try:
  pot_dir  = os.environ['POTDIR']
except:
  print 'POTDIR not defined in os environment.'
  raise

nqx = 24
nqy = 24
nqz = 24

n_sup_x = 8
n_sup_y = 8
n_sup_z = 8

pot_file  = os.path.join(pot_dir, 'PotBH.xml')
pot       = Potential('IP EAM_ErcolAd', param_filename=pot_file)
#Set paths for phonon band structures:
#paths = [np.array([0,0,-0.5]), np.array([0,0,0]), np.array([0,0,0.5]), np.array([0,0.5,0.5])]
b1        = 1.0
#Nice plot for phonon:
paths     = [np.array([0,0,0.0]), np.array([b1/2.0, b1/2.0, 0.0]), np.array([b1/2.0,b1/2.0,b1])]
paths     = calc_band_paths(paths)
unit_cell = BodyCenteredCubic(directions = [[1,0,0], [0,1,0],[0,0,1]],
                              size = (1,1,1), symbol='Fe', pbc=(1,1,1),
                              latticeconstant = 2.83)

symbols          = unit_cell.get_chemical_symbols()
cell             = unit_cell.get_cell()
scaled_positions = unit_cell.get_scaled_positions()

THZ_to_mev = 4.135665538536
unit_cell  = PhonopyAtoms(symbols=symbols, cell=cell, scaled_positions=scaled_positions)
phonon     = Phonopy(unit_cell, [[n_sup_x,0,0],[0,n_sup_y,0],[0,0,n_sup_z]]) 
phonon.generate_displacements(distance=0.05)
supercells = phonon.get_supercells_with_displacements()

#We need to get the forces on these atoms...
forces = []
print 'There are', len(supercells), 'displacement patterns'
for sc in supercells:
  cell = aseAtoms(symbols=sc.get_chemical_symbols(), scaled_positions=sc.get_scaled_positions(), 
                  cell=sc.get_cell(), pbc=(1,1,1))
  cell = Atoms(cell)
  cell.set_calculator(pot)
  forces.append(cell.get_forces())

phonon.set_forces(forces)
phonon.produce_force_constants()

mesh = [nqx, nqy, nqz]
phonon.set_mesh(mesh, is_eigenvectors=True)
qpoints, weights, frequencies, eigvecs = phonon.get_mesh()

#SHOW DOS STRUCTURE
phonon.set_total_DOS(freq_min=0.0, freq_max=12.0, tetrahedron_method=False)
phonon.get_total_DOS()
phonon.write_total_DOS()
phonon.plot_total_DOS().show()

#BAND STRUCTURE
phonon.set_band_structure(paths, is_eigenvectors=False, is_band_connection=True)
phonon.get_band_structure()
ph_plot = phonon.plot_band_structure(labels=["G", "N", "P"])
ph_plot.show()
#WRITE ANIMATION MODE
phonon.write_animation()
#TEMPERATURE
phonon.set_thermal_properties(t_step=10, t_max=1000, t_min=0)
temp_plot = open('temp_plot.dat', 'w')
for t, free_energy, entropy, cv in np.array(phonon.get_thermal_properties()).T:
    print >> temp_plot, ("%12.3f " + "%15.7f" * 3) % ( t, free_energy, entropy, cv )
phonon.plot_thermal_properties().show()
#Finally write parametrs of calculation i.e. supercell size for force_constant
#and qpoint mesh

with open('subgb.json','w') as f:
  gb_dict = {} 
  gb_dict['n_sup_cell'] = [n_sup_x, n_sup_y, n_sup_z]
  gb_dict['nqs']        = [nqx, nqy, nqz]
  json.dump(gb_dict,f,indent=2)


