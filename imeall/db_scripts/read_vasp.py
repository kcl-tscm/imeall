import os
import re
from ase.io import vasp
from quippy import Atoms
import json
import shutil

class VaspIncar(object):
  """ 
  :class:`VaspIncar` to manipulate vasp input files.
  """  
  def __init__(self, jobdir='./'):
    self.calc_type = 'relaxation'
    self.incar     = 'INCAR'
    self.kpoints   = 'KPOINTS'
    self.potcar    = 'POTCAR'
    self.poscar    = 'POSCAR'
    self.job_dir   = jobdir

  def constrained_relaxation(self,current='T   T   T', target='F   F   T'):
    """ 
    :method: `constrained_relaxation` pass target string with 
    Booleans of directions to permit relaxation along.
    """
    with open(os.path.join(self.job_dir,self.poscar),'r') as f:
      poscar = f.read()
      poscar = poscar.replace(current, target)
    with open(os.path.join(self.job_dir, self.poscar),'w') as f:
      print >> f, poscar
    return

  def cont_to_pos(self):
    """
    :method:`cont_to_pos` copy CONTCAR to POSCAR. If CONTCAR is empty POSCAR will not 
    be updated and a notification will be printed to stdout.
    """
    with open(os.path.join(self.job_dir, self.contcar),'r') as f:
      contcar = f.read()
#To determine if we have a valid CONTCAR with numbers in it:
    re_number = re.compile('[0-9]')
    if re.match(re_number, contcar):
      shutil.copy(os.path.join(self.job_dir, self.contcar), os.path.join(self.job_dir, self.contcar))
    else:
      print 'CONTCAR empty not copying files across: ', self.job_dir
    return

class VaspOutcar(object):
  """ 
  :class:`VaspOutcar` to pull information of interest
  from vasp outcar files: positions, forces, total energies.
  """  
  def __init__(self):
    self.totens = []
    self.atoms  = None
    self.A      = 0.0
    self.H      = 0.0
    self.converged = False

  def per_atom_energy(self):
    """
    Given total energy and number of atoms define per atom energy.
    """
    self.per_atom_energy = [en/float(len(self.atoms)) for en in self.toten]
    return 

  def write_json_file(self):
    """
    Create a json file for this structure.
    """
    with open('subgb.json','w') as outfile:
      cell = self.atoms.get_cell()
      A = cell[0][0]*cell[1][1]
      H = cell[2][2]
#Final energy in OUTCAR FILE.
      E_cell = self.totens[-1][2]
# For legacy purposes lets always call the relaxed energy E_gb.
      gb_dict = {'E_gb': E_cell, 
                 'converged': self.converged,
                 'A' : A, 
                 'H' : H, 'n_at': len(self.atoms), 
                 'param_file':'VASP-DFT-PBE'}
      json.dump(gb_dict, outfile, indent=2)
    return
  

  def pull_atom_forces():
    force_pos_rege = re.compile("")
    for line in pos_force.split('\n'):
      x,y,z,fx,fy,fz = map(float, line.split())

# REGEX library
# pull free energies, number of iterations, 
# and whether calculation has converged.
toten_regex = re.compile(r'free energy    TOTEN  =\s+([\-0-9\.]+)', re.S)
iter_regex  = re.compile(r'Iteration\s+([0-9]+)\(\s+([0-9]+)\)')
conv_regex  = re.compile('reached required accuracy', re.S)



if __name__=='__main__':
# Read text file
  with open('OUTCAR','r') as f:
    outcar = f.read()
  gamma_surf = VaspOutcar()
# Set whether this particular Calculation has converged.
  if conv_regex.search(outcar):
    print 'Calculation Converged'
    gamma_surf.converged = True
  
# Return an Atoms object, collect total energies and the number of iterations:
  try:
    out_atoms  = vasp.read_vasp_out('OUTCAR')
    gamma_surf.atoms = out_atoms
  except IndexError:
    print 'No position information in OUTCAR using input.'
    out_atoms  = vasp.read_vasp('POSCAR')
    gamma_surf.atoms = out_atoms

  out_atoms.write('poscar.xyz')

  totens     = map(float, toten_regex.findall(outcar))
  iterations = iter_regex.findall(outcar)
  
  if len(totens)==len(iterations):
    for step, toten in zip(iterations, totens):
      gamma_surf.totens.append((int(step[0]), int(step[1]), toten))
# If we have one more iteration than total energy
# it is likely the calculation timed out so we can just pair up to the 
# final total energy with the second last iteration
  elif len(totens) == (len(iterations)-1):
    for step, toten in zip(iterations, totens):
      gamma_surf.toten.append((int(step[0]), int(step[1]), toten))
  else:
    print "Total Energies and iterations don't agree."
    print 'Totens', len(totens), 'Iters', len(iterations)
    sys.exit()
  
# Create Total Energy Plot for all jobdirs in the pattern.
  with open('toten.dat','w' ) as f:
    for thing in gamma_surf.totens:
      print >> f, '\t'.join(map(str, thing))
  gamma_surf.write_json_file()
  gamma_surf.atoms.write('relaxed_dft_struct.xyz')
