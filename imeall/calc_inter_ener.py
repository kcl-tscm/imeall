import os
import json
import numpy as np
import argparse
from ase.optimize import FIRE
from quippy import Atoms, Potential, frange

def h2_formation_energy(pot):
  """
  Given a potential calculate the H2 formation energy and
  equilibrium bond spacing.

  Args:
    pot(:quippy:class:`Potential`) potential object.

  Returns:
    float: Hydrogen molecule formation energy.
  """
  h2 = aseAtoms('H2', positions=[[0, 0, 0],[0, 0, 0.7]])
  h2 = Atoms(h2)
  h2.set_calculator(pot)
  opt = BFGS(h2)
  opt.run(fmax=0.0001)
  E_h2  = h2.get_potential_energy()
  return E_h2

def calc_egb(json_dict):
  return json_dict['E_gb']

def get_interface_bounds(ats):
  """
  Pull the top half interface and return the boundaries
  of that interface in the original coordinates.

  Args:
    ats (:ase:class:`Atoms`): Atoms object of full bi-crystal.

  Returns: 
    gb_min, gm_max, z_width, min_at
  """
  cell_midpoint = ats.get_cell()[2,2]/2.0
  #select non BCC sites are 0 otherwise [1-3] inclusive.
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
  at_min = zint.positions[:,2].min() - gb_min
  return gb_min, gb_max, z_width, at_min

def apply_strain(ats, mode, st_num):
  """Apply a deformation mode to :class:`ase.Atoms` object. 

  Args:
    mode(str): Options are shear, stretch, hydrostatic.
    st_num(float): Strain applied as a percentage.

  Returns:
   :class:`ase.Atoms` 

  """

  e1 = np.array([1,0,0])
  e2 = np.array([0,1,0])
  e3 = np.array([0,0,1])
  cell = ats.get_cell()
  if mode == 'hydrostatic':
    strain_tensor = np.eye(3) + st_num*np.eye(3)
    cell = cell*strain_tensor
    ats.set_cell(cell, scale_atoms=True)
    print 'Hydrostatic strain', st_num
    print 'strain tensor', strain_tensor
  elif mode == 'stretch':
    strain_tensor = np.tensordot(e3, e3, axes=0)
    strain_tensor = np.eye(3) + st_num*strain_tensor
    cell = strain_tensor*cell
    print 'Stretch strain'
    print 'Cell:', cell
    ats.set_cell(cell, scale_atoms=True)
  elif mode == 'shear':
    strain_tensor = np.tensordot(e1, e2, axes=0)
    strain_tensor = np.eye(3) + st_num*strain_tensor
    cell = strain_tensor.dot(cell)
    print 'Shear Strain', strain_tensor
    print 'Cell:', cell
    ats.set_cell(cell, scale_atoms=True)
  else:
    print 'No strain applied.'
  return ats

alat = 2.83
#bulk_sites = [tetrahedral site and octahedral]
bulk_sites = map(lambda x: alat*x, [np.array([0.25, 0.0, 0.5]), np.array([0,0,0.5])])
POT_DIR = os.path.join(app.root_path, 'potentials')
eam_pot = os.path.join(POT_DIR, 'PotBH.xml')
r_scale = 1.00894848312
pot = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)


if __name__=='__main__':
  parser = argparse.ArgumentParser() 
  parser.add_argument('-m','--modes', nargs='+', help='Type of strain mode: shear, stretch, or hydrostatic.', default=['shear', 'stretch'])
  parser.add_argument('-r','--rescale', help='interstitial site needs to be rescaled into lattice', action='store_true')
  parser.add_argument('-n','--nums', nargs='+', help='Type of strain mode: shear, stretch, or hydrostatic.', default=[-0.01, -0.005, 0.0, 0.005, 0.01], type=float)
  parser.add_argument('-f','--force_tol', help='Force Tolerance', default=0.05, type=float)
  args = parser.parse_args()
  with open('unique_h_sites.json','r') as f:
    h_sites = json.load(f)
  print 'There are ', len(h_sites), 'H interstitials'
  ats = Atoms('output.xyz')
  gb_min, gb_max, z_width, at_min = get_interface_bounds(ats)
  with open('subgb.json','r') as f:
    subgb_dict = json.load(f)

  force_tol = args.force_tol
  E_h2 = -4.73831215121
  #E_gb = subgb_dict['E_gb']
  all_h_ats = ats.copy()
  for mode in args.modes:
    for num in args.nums:
      g = open('h_site_ener_{}_{}.txt'.format(mode, str(num)), 'w')
      s_ats = ats.copy()
      s_ats = apply_strain(s_ats, mode, num)
      s_ats.set_calculator(pot)
      E_gb = s_ats.get_potential_energy()
      for h_site in h_sites:
        h_ats = s_ats.copy()
        h_site_tmp = list(h_site)
        #-1.0 to subtract vacuum added in calc_eseg.py at_min to account for diff between gb_min and lowest atom.
        h_site_tmp[2] += gb_min - 1.0 + at_min
        h_ats.add_atoms(h_site_tmp, 1)
        h_ats.set_calculator(pot)
        opt = FIRE(h_ats)
        opt.run(fmax = force_tol)
        E_gbh = h_ats.get_potential_energy()
        h_at_rel = filter(lambda x:x.number == 1, h_ats)
        all_h_ats.add_atoms(h_at_rel[0].position, 1)
        #print position of relaxed h atom and the interstitial energies.
        print >> g, h_at_rel[0].position, E_gbh, E_gbh - E_gb - 0.5*E_h2
        g.flush()
      g.close()
