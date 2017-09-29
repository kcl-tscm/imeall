import os
import sys
import argparse
import subprocess
import numpy as np
from quippy import Atoms, set_fortran_indexing, Potential
from ase.constraints import UnitCellFilter, StrainFilter

set_fortran_indexing(True)

class ElasticDipole(object):
  def __init__(self):
    """Contains methods for computing
    the elastic dipole tensor of a point defect in a material.
    """
    self.strain_tensor = [-0.01, -0.005, 0.0, 0.005, 0.01]

#  def defect_force_method(self, ats, defect, rcut=None):
#    """
#    Requires an interatomic potential that can calculate
#    the distinct contribution of a force from a particular atom,
#    i.e. EAM or Tightbinding model.
#
#    Args:
#      ats (:py:class:`Atoms`): Atoms object of structure.
#      defect(:py:class:`Atom`): Specifies the defect atom.
#      rcut (float): Cutoff radius to consider forces.
#
#    Todo:
#      Implement this method
#    """
#    alpha_ij = self.calc_interatomic_force(defect_atom.index, ats)
#    return alpha_ij
    
  def relax_defect_cell(self, ats, output_name='defect_cell_relaxed.xyz', force_tol=0.0001, relax_cell=False):
    """
    Accepts atoms object with an attached calculator.Minimize forces.

    Args:
      ats (:obj:`Atoms`): Atoms with defect to relax.
      output_name (str, optional): Filename to print relaxed atoms structure to.
      force_tol (float, optional): Force tolerance to stop relaxation.
      relax_cell (bool, optional): Relax lattice parameters.

    Returns:
      :class:`Atoms`: A relaxed Atoms object.
    """
    from ase.optimize import FIRE
    if relax_cell:
      strain_mask  = [1,1,1,0,0,0]
      ucf = UnitCellFilter(ats, strain_mask)
      opt = FIRE(ucf)
    else:
      strain_mask  = [0,0,0,0,0,0]
      ucf = UnitCellFilter(ats, strain_mask)
      opt = FIRE(ucf)
    opt.run(fmax = force_tol)
    ats.write(output_name)
    return ats
  
  def compute_vacancy_dipole(self, defect, ats, pot=None, forces=np.array([])):
    """Compute dipole tensor from induced forces.

    Args: 
      defect (:py:class:`Atom`): Atom object of defect atom.
      ats (:py:class:`Atoms`): If forces are absent the :py:class:Atoms object must have defect present.
      pot (:py:class:`Potential`, optional): Potential for calculating interatomic forces if required.
      forces(:py:class:`numpy.array`, optional): numpy array of forces if already available.

    Returns: 
      3x3 numpy array of G the dipole tensor.
    """
    if not forces.any():
      print len(ats) 
      ats.remove_atoms([defect.index+1])
      ats.set_calculator(pot)
      forces = ats.get_forces()
      ats.write('relaxed_cell_removed_defect.xyz')
    else:
      #force array has been passed (e.g. read from OUTCAR)
      pass
    alpha_ij = np.zeros([3,3])
    for defect_force, at in zip(forces, ats):
      if np.linalg.norm(defect.position-at.position) <= 'Inf': 
        alpha_ij[0, 0] += defect_force[0]*(defect.position[0] - at.position[0])
        alpha_ij[0, 1] += defect_force[0]*(defect.position[1] - at.position[1])
        alpha_ij[0, 2] += defect_force[0]*(defect.position[2] - at.position[2])
        alpha_ij[1, 0] += defect_force[1]*(defect.position[0] - at.position[0])
        alpha_ij[1, 1] += defect_force[1]*(defect.position[1] - at.position[1])
        alpha_ij[1, 2] += defect_force[1]*(defect.position[2] - at.position[2])
        alpha_ij[2, 0] += defect_force[2]*(defect.position[0] - at.position[0])
        alpha_ij[2, 1] += defect_force[2]*(defect.position[1] - at.position[1])
        alpha_ij[2, 2] += defect_force[2]*(defect.position[2] - at.position[2])
    return alpha_ij


def find_h_atom(ats):
  """Finds the hydrogen atom in .xyz file and returns its :py:class:`Atom` object.
  
  Args:
    ats(:py:class:`Atoms`)

  Returns: 
    :py:class:`Atom` object.
  """
  h_list = [at for at in ats if at.number==1]
  if len(h_list) > 1:
    sys.exit('too many hydrogens in unit cell.')
  elif len(h_list) == 0:
    sys.exit('No hydrogen in unit cell.')
  else:
    defect = h_list[0]
  return defect

def calc_elast_dipole_eam(input_file, force_tol, relax_cell):
  """
  Calculate elastic dipole using an Embedded Atom Potential.

  Args:
    input_file (str): Name of .xyz file contains unitcell with defect.
    force_tol (float): Force tolerance to stop relaxation.
    relax_cell (bool): Relax lattice vectors.

  Returns:
    Elastic Dipole Tensor 3x3 numpy array.
  """
  try:
    POT_DIR = os.environ['POTDIR']
  except KeyError:
    sys.exit("PLEASE SET export POTDIR='path/to/potfiles/'")

  elastic = ElasticDipole()  

  eam_pot = os.path.join(POT_DIR, 'PotBH.xml')
  pot = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(1.00894848312), param_filename=eam_pot)

  ats = Atoms(input_file)
  ats.set_calculator(pot)

  init_vol = ats.get_volume()
  print 'Initial Vol.', init_vol
  elastic.relax_defect_cell(ats, force_tol=force_tol, relax_cell=relax_cell)
  final_vol = ats.get_volume()
  print 'Final Vol.', final_vol
  print 'Volume Diff.', final_vol-init_vol
  ats = Atoms('defect_cell_relaxed.xyz')
  defect = find_h_atom(ats)
  print 'Defect index', defect.index, 'Position', defect.position, 'Type: ', defect.number
  ats.set_calculator(pot)
  return elastic.compute_vacancy_dipole(defect, ats.copy(), pot)

def calc_elast_dipole_dft(input_file):
  """Reads OUTCAR file in the same directory with one shot forces 
  induced by removal of defect. Reads defect position 
  from .xyz file (which contains the defect) defined by `input_file` 
  calculates and returns the elastic dipole tensor of the defect.

  Args:
    input_file(str): name of input .xyz file containing defect cell.

  Returns:
    Elastic Dipole Tensor 3x3 numpy array.
  """
  import ase.io.vasp as vs
  elastic = ElasticDipole()
  ats_def = Atoms(input_file)
  defect = find_h_atom(ats_def)
  ats_pos = vs.read_vasp()
  ats = vs.read_vasp_out()
  ats = Atoms(ats)
  print 'Defect index', defect.index, 'Position', defect.position, 'Type: ', defect.number
  ats.get_forces()
  ats.write('force.xyz')
  return elastic.compute_vacancy_dipole(defect, ats_pos.copy(), forces=ats.get_forces())

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-ct', '--calc_type', help="Can use 'EAM' or 'DFT' potential. If DFT then script requires OUTCAR with single shot \
                                                  forces for defectless cell and xyz input for cell with defect.", required=True)
  parser.add_argument('-rc', '--relax_cell', help="Relax defect super cell.", action='store_true')
  parser.add_argument('-f', '--force_tol', help="force_tolerance.", type=float, default=0.0001)
  parser.add_argument('-i', '--input_file', help="Structure file of system with defect. (.xyz format)", required=True)
  args = parser.parse_args()

  if args.calc_type == 'EAM':
    print calc_elast_dipole_eam(args.input_file, args.force_tol, args.relax_cell).round(3)
  elif args.calc_type == 'DFT':
    print calc_elast_dipole_dft(args.input_file)
  else:
    sys.exit('Invalid potential type chosen. Please select DFT or EAM.')
