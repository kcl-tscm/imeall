import os
import re
from ase.io import vasp
from quippy import Atoms
import json

class VaspOutcar(object):
  """ 
  :class:`VaspOutcar` to pull information of interest
  from vasp calculation.
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
  out_atoms  = vasp.read_vasp_out('OUTCAR')
  gamma_surf.atoms = out_atoms
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
  
# Create Total Energy Plot
  with open('toten.dat','w' ) as f:
    for thing in gamma_surf.totens:
      print >> f, '\t'.join(map(str, thing))
  gamma_surf.write_json_file()
