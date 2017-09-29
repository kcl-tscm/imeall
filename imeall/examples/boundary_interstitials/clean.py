import os

files_to_remove = ['POSCAR',
'h_energetics.xyz',
'h_energetics.xyz.idx',
'h_site_ener_stretch_0.0.txt',
'hydrogenated_grain.xyz',
'imeall.log',
'interface.xyz',
'interface.xyz.idx',
'output.xyz.idx',
'unique_h_sites.json',
'unique_lattice_sites.json']

for f in files_to_remove:
  try:
    os.remove(f)
  except OSError:
    pass
