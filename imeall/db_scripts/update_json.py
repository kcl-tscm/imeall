import os
import sys
import json
from quippy import Atoms


jobdirs = []
for thing in os.listdir('./'):
  if os.path.isdir(thing) and thing[:3]=='111':
    jobdirs.append(thing)

for job in jobdirs[2:]:
  json_path = os.path.join(job, 'gb.json')
  new_json = {}
  with open(json_path,'r') as json_old:
    old_json = json.load(json_old) 
    new_json['zplanes'] = old_json['zplanes']
    new_json['orientation_axis'] = old_json['orientation axis']
    new_json['boundary_plane']   = old_json['boundary plane']
    new_json['coincident_sites'] = old_json['coincident sites']
    new_json['angle'] = old_json['angle']
    new_json['gbid']  = old_json['gbid']
    new_json['n_at']  = old_json['n_unit_cell']
    new_json['type']  = 'symmetric tilt boundary'
    at = Atoms('{0}.xyz'.format(os.path.join(job, (old_json['gbid']))))    
    cell = at.get_cell()
    A    = cell[0,0]*cell[1,1]
    new_json['A']  = A

  with open(json_path,'w') as json_new_file:
    json.dump(new_json, json_new_file, indent=2)


