import os
import numpy as np
from ase import Atoms as aseAtoms
from quippy   import Atoms, Potential, Atoms
from phonopy  import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from ase.lattice.cubic       import BodyCenteredCubic


try:
  pot_dir  = os.environ['POTDIR']
except:
  print 'POTDIR not defined in os environment.'
  raise

pot_file = os.path.join(pot_dir, 'PotBH.xml')
pot     = Potential('IP EAM_ErcolAd', param_filename=pot_file)

unit_cell = BodyCenteredCubic(directions = [[1,0,0], [0,1,0],[0,0,1]],
                              size = (1,1,1), symbol='Fe', pbc=(1,1,1),
                              latticeconstant = 2.83)
symbols          = unit_cell.get_chemical_symbols()
cell             = unit_cell.get_cell()
scaled_positions = unit_cell.get_scaled_positions()


n_sup_cell = 8
THZ_to_mev = 4.135665538536
unit_cell  = PhonopyAtoms(symbols=symbols, cell=cell, scaled_positions=scaled_positions)
phonon     = Phonopy(unit_cell, [[n_sup_cell,0,0],[0,n_sup_cell,0],[0,0,n_sup_cell]]) 
phonon.generate_displacements(distance=0.1)
supercells = phonon.get_supercells_with_displacements()
#We need to get the forces on these atoms...
forces = []
for sc in supercells:
  cell = aseAtoms(symbols=sc.get_chemical_symbols(), scaled_positions=sc.get_scaled_positions(), 
                  cell=sc.get_cell(), pbc=(1,1,1))
  cell = Atoms(cell)
  cell.set_calculator(pot)
  forces.append(cell.get_forces())
phonon.set_forces(forces)
print phonon.produce_force_constants()
mesh = [20, 20, 20]
phonon.set_mesh(mesh, is_eigenvectors=True)
qpoints, weights, frequencies, eigvecs = phonon.get_mesh()
#SHOW DOS STRUCTURE
phonon.set_total_DOS(freq_min=0.0, freq_max=60.0)
phonon.get_total_DOS()
phonon.write_total_DOS()
phonon.plot_total_DOS().show()
phonon.write_animation(band_index=1)
#TEMPERATURE
phonon.set_thermal_properties(t_step=10, t_max=1000, t_min=0)
temp_plot = open('temp_plot.dat', 'w')
for t, free_energy, entropy, cv in np.array(phonon.get_thermal_properties()).T:
    print >> temp_plot, ("%12.3f " + "%15.7f" * 3) % ( t, free_energy, entropy, cv )
phonon.plot_thermal_properties().show()
