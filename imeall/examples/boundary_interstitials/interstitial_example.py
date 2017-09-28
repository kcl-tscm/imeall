from imeall.calc_eseg import  gen_interface
from imeall.calc_inter_dist import decorate_interface
from imeall.calc_inter_ener import gen_interface_bounds
from imeall import app


#can calculate dissolution energy for different strain modes.
mode = 'stretch' 
num = 0.0
E_h2 = -4.73831215121
force_tol = 0.1

POT_DIR = os.path.join(app.root_path, 'potentials')
eam_pot = os.path.join(POT_DIR, 'PotBH.xml')
r_scale = 1.00894848312
pot = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)

if __name__=='__main__':
  gen_interface()
  decorate_interface()
  with open('unique_h_sites.json','r') as f:
    h_sites = json.load(f)
  g = open('h_site_ener_{}_{}.txt'.format(mode, str(num)), 'w')
  ats = Atoms('output.xyz')
  gb_min, gb_max, z_width, at_min = get_interface_bounds(ats)
  with open('subgb.json','r') as f:
    subgb_dict = json.load(f)
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

