import os
import json
import numpy as np
from ase.optimize import FIRE
from quippy import Atoms, Potential, frange


#calc_inter_diss_ener.py calculates intersitial dissolution energy
def h2_formation_energy(pot):
  h2 = aseAtoms('H2', positions=[[0, 0, 0],[0, 0, 0.7]])
  h2 = Atoms(h2)
  h2.set_calculator(pot)
  opt = BFGS(h2)
  opt.run(fmax=0.0001)
  E_h2  = h2.get_potential_energy()

def calc_egb(json_dict):
  return json_dict['E_gb']

def get_interface_bounds(ats):
  """
  Pull the top half interface and return the boundaries
  of that interface in the original coordinates.
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
  return gb_min, gb_max, z_width

#calculated for Ramasubramaniam

#Calculates dissolution energies given a list of h_sites.
#bulk_sites = [tetrahedral site and octahedral]
alat = 2.83
bulk_sites = map(lambda x: alat*x, [np.array([0.25, 0.0, 0.5]), np.array([0,0,0.5])])

POT_DIR = os.environ['POTDIR']
eam_pot = os.path.join(POT_DIR, 'PotBH.xml')
r_scale = 1.00894848312
pot = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)

with open('unique_h_sites.json','r') as f:
  h_sites = json.load(f)

ats = Atoms('output.xyz')
gb_min, gb_max, z_width = get_interface_bounds(ats)
g = open('h_site_ener.txt','w')
#CONSTANTS
with open('subgb.json','r') as f:
  subgb_dict = json.load(f)

force_tol = 0.01
E_h2 = -4.73831215121
#E_h2_dft = -4.52
E_gb = subgb_dict['E_gb']
#Test
all_h_ats = ats.copy()
for h_site in h_sites:
  h_site[2] += float(gb_min)-1.00000 #rescale into cell remove vacuum of 1 A.
  h_ats = ats.copy()
  if (gb_min + 1.0*z_width) <= h_site[2] <= (gb_max-1.0*z_width): #test tight to interface
    all_h_ats.add_atoms(h_site,1)
    h_ats.add_atoms(h_site,1)
    h_ats.set_calculator(pot)
    opt = FIRE(h_ats)
    opt.run(fmax = force_tol)
    E_gbh = h_ats.get_potential_energy()
    h_at_rel = filter(lambda x:x.number == 1, h_ats)
    print h_at_rel
    print >> g, h_at_rel[0].position, E_gbh, E_gbh - E_gb - 0.5*E_h2
    g.flush()
  else:
    pass

for i in frange(len(h_ats)):
  all_h_ats.id[i] = i
all_h_ats.write('full_hydrogenated_grain.xyz')
g.close()
