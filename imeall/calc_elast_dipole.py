import os
import argparse
import subprocess
from ase.constraints import UnitCellFilter, StrainFilter
from quippy import Atoms, set_fortran_indexing

set_fortran_indexing(False)

class ElasticDipole(object):
  def __init__(self,):
    self.strains=[−0.01, −0.005, 0, 0.005, 0.01]

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
      ats: `Atoms` object of system.
      defect: `Atom` object specifying the defect atom.
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
    strain_mask  = [1,1,1,0,0,0]
    if relax_cell:
      ucf = UnitCellFilter(ats, strain_mask)
      opt = FIRE(ucf)
    else:
      opt = FIRE(ats)
    opt.run(fmax = force_tol)
    ats.write(name)
    return ats
  
  def compute_vacancy_dipole(defect, pot, input_cell='defect_cell_relaxed.xyz'):
    ats = Atoms(input_cell)
    ats.set_calculator(pot)
    del(ats[defect.index])
    forces = ats.get_atomic_forces()
    #calculat dipole tensor
    alpha_ij = np.zeros(3,3)
    for at in ats:
      for k in range(3):
        if at.index != defect.index:
          self.calc_interatomic_force(at, defect_atom)
          alpha_ij[0, 0] = defect_force[0]*(at.position[0]-defect.position[0])
          alpha_ij[0, 1] = defect_force[0]*(at.position[1]-defect.position[1])
          alpha_ij[0, 2] = defect_force[0]*(at.position[2]-defect.position[2])
          alpha_ij[1, 0] = defect_force[1]*(at.position[0]-defect.position[0])
          alpha_ij[1, 1] = defect_force[1]*(at.position[1]-defect.position[1])
          alpha_ij[1, 2] = defect_force[1]*(at.position[2]-defect.position[2])
          alpha_ij[2, 0] = defect_force[2]*(at.position[0]-defect.position[0])
          alpha_ij[2, 1] = defect_force[2]*(at.position[1]-defect.position[1])
          alpha_ij[2, 2] = defect_force[2]*(at.position[2]-defect.position[2])
    return alpha_ij

if __name__='__main__':
  parser = argparse.ArgumentParser()
  parser.add_argumet('-i', '--input_file', help="Structure file of system with defect. (xyz)", required=True)
  parser.add_argumet('-d', '--defect_index', help="Atom index (base 0) of the defect atom.", require=True)

  try:
    POT_DIR     = os.environ['POTDIR']
  except:
    sys.exit("PLEASE SET export POTDIR='path/to/potfiles/'")
  elastic = ElasticDipole()  
  eam_pot = os.path.join(POT_DIR, 'PotBH.xml')
  pot = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(1.00894848312), param_filename=eam_pot)

# relax defect cell
  ats = Atoms(args.input_file)
  defect_index = 0
  defect = ats[defect_index]
  elastic.relax_defect_cell(ats)
# compute dipole tensor
  elastic.compute_vacancy_dipole(defect)
