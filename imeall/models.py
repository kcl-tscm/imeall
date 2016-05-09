import os
import sys
from flask  import Flask, request, session, g, redirect
from  flask import    url_for, abort, render_template, flash
from   quippy import Atoms
from   imeall import app
import imeall.slabmaker.slabmaker as slabmaker
import json

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

  def update_json(self, dirname):
    ''' This function was originally written to update all keys in the
    json dictionaries in the grain boundary directories.
    The pattern is quite general and can be adapted to just add
    new keys delete old keys consider it a dictionary migration
    routine.'''
    os.path.join(dirname,'gb.json')
    new_json
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

