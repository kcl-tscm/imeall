import json

from quippy import Atoms
from imeall.calc_inter_ener import get_interface_bounds

with open('unique_h_sites.json','r') as f:
    h_sites = json.load(f)
print 'There are ', len(h_sites), 'H interstitials'

ats = Atoms('output.xyz')
gb_min, gb_max, z_width, at_min = get_interface_bounds(ats)
h_ats = ats.copy()

for h_site in h_sites:
    h_site_tmp = list(h_site)
    #remove vacuum restore min position
    h_site_tmp[2] += gb_min - 1.0 + at_min
    h_ats.add_atoms(h_site_tmp, 1)

h_int_ats = Atoms('interface.xyz')
for h_site in h_sites:
    h_site_tmp = list(h_site)
    h_int_ats.add_atoms(h_site_tmp, 1)

for i in range(0,len(h_ats)):
    h_ats.id[i] = i

for i in range(0, len(h_int_ats)):
    h_int_ats.id[i] = i

h_ats.write('decorated.xyz')
h_int_ats.write('int_decorated.xyz')
