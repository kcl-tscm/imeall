import os
import re
import sys
import json
import numpy as np
from   quippy import Atoms
try:
  import imeall.slabmaker.slabmaker as slabmaker
  from  flask  import Flask, request, session, g, redirect
  from  flask import    url_for, abort, render_template, flash
except:
  print 'no flask'

# Currently Our models are stored by hand
# and then we handle the interactions with 
# directory structure manually. Which should allow us to serve
# files.

class Job(object):
  def __init__(self):
    self.pbs_file = ''
    self.job_dir  = ''
    self.job_id   = ''

  def sub_pbs(self, job_dir, exclude='DFT', suffix='v6bxv2z', regex=None):
    """ 
    Given an explicit suffix, or a regex this routine recurses through
    the directory structure and submits any pbs files that 
    match the suffix or regex pattern. Exclude keeps track of 
    directories that (we mightn't want for instance DFT on Ada
    or EAM on Mira.)
    Useful submission patterns include:
    REGEX:
      submit all pbs files with a rigid body translation 
      and any atom deletion criterion hydrogen concentration etc.: 
        tv[.0-9]+bxv[.0-9]+_.*?
      submit all sub-pbs files with a rigid body translation
      and different atom _deletion criterion: 
        tv[.0-9]+bxv[.0-9]+_d[.0-9]+
      all translations with a specific deletion criterion
      in this case 2.3 A:
        tv[.0-9]+bxv[.0-9]+_d2.3
      etc.
    SUFFIX:
      submit all super cells: 
        v6bxv2z
    """
    lst = os.listdir(job_dir)
    for dir in lst:
      dir = os.path.join(job_dir, dir)
      if regex == None:
        if os.path.isdir(dir) and dir != 'DFT':
          print dir
          self.sub_pbs(dir, suffix=suffix, regex=regex)
        elif dir.split('_')[-1] == suffix:
          pbs_dir = os.path.join(sub_dir, dir)
          os.system("cd {0}; qsub fe{1}.pbs".format(pbs_dir, job_dir+'_'+suffix))
        else:
          pass
      else:
        if os.path.isdir(dir) and dir != 'DFT':
          self.sub_pbs(dir, suffix=suffix, regex=regex)
        elif regex.match(dir):
          try:
            dir  = '/'.join(dir.split('/')[:-1])
            name = dir.split('/')[-1]
            os.system("cd {0}; qsub fe{1}.pbs".format(dir, name))
          except:
            print 'Job Submit Failed'
        else:
          pass

class GBMaintenance(object):
  """
  Collection of maintenance routines for the GB database.
  Possible usages: regenerate all the csl lattices in the database
  or a subdirectory of the database, take a new grain boundary profile 
  picture for multiple directories, update the gb json information if a
  new grain boundary property is desired. This really is 
  turning into facebook for grain boundaries!
  """
  def __init__(self):
    self.materials = ['alphaFe']

  def retake_pic(self,fname, translate=False,toggle=False, confirm=True):
    """ 
    Take grain boundary profile pic in directory
    requires gb directory with gbid.xyz file in it.
    set confirm = False to not prompt for overwrite.
    """
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

  def remove_eo_files(self, path):
    """
    In case the rsync brings across a bunch of log files
    we can get rid of those.
    """
    eo_regex = re.compile(r'[eo][0-9]+')
    lst = os.listdir(path)
    for filename in lst:
      filename = os.path.join(path, filename)
      if os.path.isdir(filename):
        self.remove_eo_files(filename)
      elif eo_regex.match(filename.split('.')[-1]):
        print filename
        os.remove(filename)
      else:
        pass

  def remove_xyz(self, path, dryrun=True, remove_type='traj'):
    """ 
    Removes xyz files according to a particular regex
    if test is true it only prints the files to be removed.
    """
    if remove_type == 'traj':
      xyz_regex = re.compile(r'.*?traj.*?xyz')
    elif remove_type == 'rbt': 
      xyz_regex = re.compile(r'.*?tv[.0-9]+bxv[.0-9]+.xyz')
    elif remove_type == 'structs': 
      xyz_regex = re.compile(r'.*?tv[.0-9]+bxv[.0-9]+_d[.0-9]+z.xyz')
    elif remove_type == 'all': 
      xyz_regex = re.compile(r'.*?xyz')
    else:
      sys.exit("No type chosen")
    lst = os.listdir(path)
    for filename in lst:
      filename = os.path.join(path, filename)
      if os.path.isdir(filename):
        self.remove_xyz(filename, dryrun = dryrun, remove_type=remove_type)
      elif xyz_regex.match(filename):
        if dryrun == True:
      #could add a function here which opens the subgb.json file 
      #checks if it is converged if not do not delete:
          dir_name = os.path.dirname(filename)
          subgb_file = os.path.join(dir_name,'subgb.json')
          if os.path.exists(subgb_file):
            with open(subgb_file,'r') as f:
              subgb_dict = json.load(f)
            if 'converged' in subgb_dict:
              print 'Converged: ', subgb_dict['converged']
          print filename, os.path.getsize(filename)
        elif dryrun == False:
          dir_name = os.path.dirname(filename)
          subgb_file = os.path.join(dir_name,'subgb.json')
          if os.path.exists(subgb_file):
            with open(subgb_file,'r') as f:
              subgb_dict = json.load(f)
            if 'converged' in subgb_dict:
              if subgb_dict['converged']:
                print 'Removing', filename
                os.remove(filename)
              else:
                print 'Not converged leaving xyz', filename
                pass
          else:
            os.remove(filename)
      else:
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
    """ 
    This function was originally written to update all keys in the
    json dictionaries in the grain boundary directories.
    The pattern is quite general and can be adapted to just add
    new keys, delete old keys, consider it a dictionary migration
    routine.
    """
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

  def fix_json(self, path):
    """
    Once my json files had two {}{} dictionaries written to them
    this parser opened all the subgb files, 
    and selected the dictionary I actually wanted.
    """
    lst = os.listdir(path)
    for filename in lst:
      new_path = os.path.join(path, filename)
      if os.path.isdir(new_path):
	      fix_json(new_path)
      elif new_path[-10:] == 'subgb.json':
	      try: 
	        with open(new_path,'r') as f:
	          j_file = json.load(f)
	      except ValueError:
	        print 'Value Error', new_path
	        with open(new_path,'r') as f:
	          j_file = f.read()
	        with open(new_path,'w') as f:
	          print >> f, j_file.split('}')[1]+'}'
	      try:
	        with open(new_path,'r') as f:
	          j_file = json.load(f)
	        print 'j_file fixed'
	      except:
	        print new_path, 'Still Broken'
      else:
        pass

  def update_json(self, filename):
    """ 
    This function was originally written to update all keys in the
    json dictionaries in the grain boundary directories.
    The pattern is quite general and can be adapted to just add
    new keys delete old keys consider it a dictionary migration
    routine.
    """
    new_json = {}
    with open(filename,'r') as json_old:
      old_json = json.load(json_old)
      new_json['zplanes'] = old_json['zplanes']
      new_json['orientation_axis'] = old_json['orientation axis']
      new_json['boundary_plane']   = old_json['boundary plane']
      new_json['coincident_sites'] = old_json['coincident sites']
      new_json['angle'] = old_json['angle']
      new_json['gbid']  = old_json['gbid']
      new_json['n_at']  = old_json['n_unit_cell']
      new_json['type']  = 'symmetric tilt boundary'
      dir_path = os.path.join('/'.join((filename.split('/'))[:-1]), old_json['gbid'])
      at = Atoms('{0}.xyz'.format(dir_path, old_json['gbid']))
      cell = at.get_cell()
      A    = cell[0,0]*cell[1,1]
      new_json['A']  = A
      json_path = filename
    with open(json_path,'w') as json_new_file:
      json.dump(new_json, json_new_file, indent=2)


class GBAnalysis():
  def __init__(self):
    pass

  def find_gb_json(self, path, j_list, filetype):
    """ 
    Returns list of directories containing grain json files
    and the filename of the json files.
    """
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

  def extract_energies(self, or_axis='001'):
#   pull GB formation energies in two stage recursive process
#   go into a grain boundary directory, recurse down through
#   grain boundary to find gb_energies pull them out and plot them
#   returns dictionary []
#   the database should only contain unique grain boundaries
#   so no key should be overwritten.
    gb_files = []
    self.find_gb_json('./imeall/grain_boundaries/alphaFe/{0}/'.format(or_axis), gb_files, 'gb.json')
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
          gb_ener = 16.02*((sub_dict['E_gb']-(-4.2731*sub_dict['n_at']))/(2*sub_dict['A']))
          #gb_ener = 16.02*((sub_dict['E_gb_init']-(-4.2731*sub_dict['n_at']))/(2*sub_dict['A']))
          grain_energy_dict['energies'].append(gb_ener)
        except:
          print 'Couldnt Extract GB energy'
      grain_energies.append(grain_energy_dict)
    return grain_energies

if __name__ == '__main__':
  print 'Should Write Some Error/Consistency Checking Code here.'
