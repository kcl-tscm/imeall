import os
import sys
from flask  import Flask, request, session, g, redirect
from  flask import    url_for, abort, render_template, flash
from   quippy import Atoms
from   imeall import app
import imeall.slabmaker.slabmaker as slabmaker
import json
import numpy as np

# Currently Our models are stored by hand
# and then we handle the interactions with 
# directory structure manually. Which should allow us to serve
# files.

class GBMaintenance(object):
  '''
  Collection of maintenance routines for the GB database.
   Possible usages: regenerate all the csl lattices in the database
   or a subdirectory of the database, take a new grain boundary profile 
   picture for multiple directories, update the gb json information if a
   new grain boundary property is desired. This really is 
   turning into facebook for grain boundaries!'''
  def __init__(self):
    self.materials = ['alphaFe']

  def retake_pic(self,fname, translate=False,toggle=False, confirm=True):
    ''' retake grain boundary profile pic in directory
      requires gb directory with gbid.xyz file in it.
      set confirm = False to not prompt for overwrite'''
    if confirm:
      var = 'n'
      var = raw_input('Retake photo (y/n)?')
    else:
      var = 'y'

    if var =='y':
      fname = os.path.join(fname,fname)
      slabmaker.take_pic(fname, translate=translate, toggle=toggle)
      print 'retaking photo'
    elif var =='n':
      pass

  def add_key_to_dict(self, dirname):
    os.path.join(dirname, 'subgb.json')
    new_json = {}
    with open(json_path,'r') as json_old:
      old_json = json.load(json_old)
    for key in old_json.keys():
      new_json[key] = old_json[key]
    at = Atoms('{0}.xyz'.format(os.path.join(dirname, )))
    cell = at.get_cell()
    A    = cell[0,0]*cell[1,1]
    new_json['A']    = A
    new_json['n_at'] = len(at) 

  def update_json(self, dirname):
    ''' This function was originally written to update all keys in the
    json dictionaries in the grain boundary directories.
    The pattern is quite general and can be adapted to just add
    new keys delete old keys consider it a dictionary migration
    routine.'''
    os.path.join(dirname,'gb.json')
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


class GBAnalysis():
  def __init__(self):
    pass

  def find_gb_json(self, path, j_list, filetype):
    ''' 
    Returns list of directories containing grain json files
    and the filename of the json files.
    '''
    try:
      lst = os.listdir(path)
    except:
      pass
    for filename in lst:
      filename = os.path.join(path,filename)
      if os.path.isdir(filename) and filename.split('/')[-1] != 'DFT':
        self.find_gb_json(filename, j_list, filetype)
      elif filename.split('/')[-1] == filetype:
        j_list.append([path,filename])
      else:
        pass

  def extract_energies(self):
#   pull GB formation energies in two stage recursive process
#   go into a grain boundary directory, recurse down through
#   grain boundary to find gb_energies pull them out and plot them
#   returns dictionary []
#   the database should only contain unique grain boundaries
#   so no key should be overwritten.
    gb_files = []
    self.find_gb_json('./grain_boundaries/alphaFe/110/', gb_files, 'gb.json')
    grain_energies = []
    for grain in  gb_files:
      grain_energy_dict = {}
      path = grain[0]
      with open(grain[1],'r') as f:
        j_dict = json.load(f)
      grain_energy_dict['orientation_axis'] = j_dict['orientation_axis']
      grain_energy_dict['angle']            = j_dict['angle']*(180./np.pi)
      grain_energy_dict['boundary_plane']   = j_dict['boundary_plane']
      grain_energy_dict['energies']         = []
      subgb_files = []
      self.find_gb_json(path, subgb_files, 'subgb.json')
      for subgrain in subgb_files:
        with open(subgrain[1],'r') as f:
          sub_dict = json.load(f)
        try:
          #gb_ener = 16.02*((sub_dict['E_gb']-(-4.2731*sub_dict['n_at']))/(2*sub_dict['A']))
          gb_ener = 16.02*((sub_dict['E_gb_init']-(-4.2731*sub_dict['n_at']))/(2*sub_dict['A']))
          grain_energy_dict['energies'].append(gb_ener)
        except:
          print 'Error on Atoms'
        #except KeyError:
        #  print subgrain[0], 'Key Missing'
      grain_energies.append(grain_energy_dict)
          
    return grain_energies

if __name__ == '__main__':
  analyze =  GBAnalysis()
  gb_list = analyze.extract_energies()
  print '0.0 0.0 0.0 0.0'
  for gb in sorted(gb_list, key = lambda x: x['angle']):
    print gb['angle'], ' '.join(map(str, gb['energies']))
  print '180.0 0.0 0.0 0.0'
 
