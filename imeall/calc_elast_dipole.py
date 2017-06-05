import os
import argparse
import subprocess
import numpy as np
from quippy import Atoms, set_fortran_indexing, Potential
from ase.constraints import UnitCellFilter, StrainFilter

set_fortran_indexing(False)

class ElasticDipole(object):
  def __init__(self):
    """
    Default numerical differentiation.
    """
    self.strain_tensor = [-0.01, -0.005, 0.0, 0.005, 0.01]

  def calc_interatomic_force(self, defect_index, at_index):
    """
    Should be in QUIP.
    """
    pass

  def defect_force_method(self, ats, defect, rcut=None):
    """
    Requires an interatomic potential that can calculate
    the distinct contribution of a force from a particular atom,
    i.e. EAM or Tightbinding model.
    Params:
      :ats: `Atoms` object of system.
      :defect: `Atom` object specifying the defect atom.
    """
    alpha_ij = self.calc_interatomic_force(defect_atom.index, ats)
    return alpha_ij
    
  def relax_defect_cell(self, ats, output_name='defect_cell_relaxed.xyz', force_tol=0.0001, relax_cell=False):
    """
    Accepts atoms object with an attached calculator.Minimize forces.
    Params:
      :ats: atoms
      :force_tol: force tolerance
      :relax_cell: relax lattice parameters
    """
    from ase.optimize import FIRE
    if relax_cell:
      strain_mask  = [1,1,1,0,0,0]
      ucf = UnitCellFilter(ats, strain_mask)
      opt = FIRE(ucf)
    else:
      opt = FIRE(ats)
    opt.run(fmax = force_tol)
    ats.write(output_name)
    return ats
  
  def compute_vacancy_dipole(self, defect, ats, pot=None, forces=None):
    """
    :method: Valid where there is a clear distinction between host lattice
    and the defect atom.
    """
    if forces == None:
      ats.remove_atoms(defect.index+1)
      ats.set_calculator(pot)
      forces = ats.get_forces()
      ats.write('relaxed_cell_removed_defect.xyz')
    else:
      #force array has been passed (e.g. read from OUTCAR)
      pass
    alpha_ij = np.zeros([3,3])
    for defect_force, at in zip(forces, ats):
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

if __name__=='__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-i', '--input_file', help="Structure file of system with defect. (xyz)", required=True)
  parser.add_argument('-d', '--defect_index', help="Atom index (base 0) of the defect atom.", type=int)
  parser.add_argument('-rc', '--relax_cell', help="Relax defect super cell.", action='store_true')
  parser.add_argument('-f', '--force_tol', help="force_tolerance.", type=float, default=0.0001)
  args = parser.parse_args()
  try:
    POT_DIR     = os.environ['POTDIR']
  except:
    sys.exit("PLEASE SET export POTDIR='path/to/potfiles/'")

  elastic = ElasticDipole()  
  eam_pot = os.path.join(POT_DIR, 'PotBH.xml')
  pot = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(1.00894848312), param_filename=eam_pot)

# relax defect cell
  ats = Atoms(args.input_file)
  ats.set_calculator(pot)

  init_vol = ats.get_volume()
  print 'Initial Vol.', init_vol
  elastic.relax_defect_cell(ats, force_tol=args.force_tol, relax_cell=args.relax_cell)
  final_vol = ats.get_volume()
  print 'Final Vol.', final_vol
  print 'Diff.', final_vol-init_vol
# compute dipole tensor
  ats = Atoms('defect_cell_relaxed.xyz')
  h_list = [at for at in ats if at.number==1]
  if len(h_list) > 1:
    sys.exit('too many hydrogens in unit cell')
  else:
    defect = h_list[0]
  print 'Defect index', defect.index, 'Position', defect.position, 'Type: ', defect.number
  ats.set_calculator(pot)
  print elastic.compute_vacancy_dipole(defect, ats.copy(), pot)
