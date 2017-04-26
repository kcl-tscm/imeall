import os
import re
import sys
import glob
import json
import logging
import argparse
import numpy as np
import slabmaker.slabmaker as slabmaker

from  quippy import Atoms
from  scipy.spatial import Voronoi, voronoi_plot_2d

try:
  from  imeall import app
except ImportError:
  print 'No app'

try:
  from flask  import Flask, request, session, g, redirect
  from flask  import url_for, abort, render_template, flash
except:
  print 'No Flask Module Available'

# Currently Our models are stored by hand
# and then we handle the interactions with 
# directory structure manually. Which should allow us to serve
# files.

class PotentialParameters(object):
  """
  :class:`PotentialParameters` contains the ground state energy per atom for bulk iron
  of the different EAM potentials used in the Imeall Database. It also contains
  the values required to rescale the relaxed lattice parameters for bulk iron predicted by
  the EAM potentials to the target DFT value.
  """
  def __init__(self):
    self.name = 'Potential Parameters'

  def gs_ener_per_atom(self): 
    eperat = {'Fe_Mendelev.xml' : -4.12243503431,
              'PotBH.xml'       : -4.01298214176,
              'iron_mish.xml'   : -4.28000356875,
              'Fe_Ackland.xml'  : -4.01298226805,
              'Fe_Dudarev.xml'  : -4.31608690638,
              'dft_vasp_pbe'    : -8.238035,
              'gp33b.xml'       : -3460.93341688
             }
    return eperat

  def eam_rscale(self): 
    rscale = {'Fe_Mendelev.xml' : 1.00894848312,
              'PotBH.xml'       : 1.00894848312,
              'iron_mish.xml'   : 1.0129007626,
              'Fe_Ackland.xml'  : 1.00894185389,
              'Fe_Dudarev.xml'  : 1.01279093417,
              'dft_vasp_pbe'    : 1.00000000000,
              'gp33b.xml'       : 1.0015226318
              }
    return rscale

  def paramfile_dict(self):
    """
    Dictionary Mapping directories to quip xml files.
    Returns: 
      {potential_directory:potential_file}
    """
    paramfile      = {'DFT':'dft_vasp_pbe',
                      'PotBH':'PotBH.xml',
                      'EAM_Ack':'Fe_Ackland.xml',
                      'EAM_Men':'Fe_Mendelev.xml',
                      'EAM_Mish':'iron_mish.xml',
                      'EAM_Dud':'Fe_Dudarev.xml',
                      'GAP':'gp33b.xml'}
    return paramfile

  def potdir_dict(self):
    """
    :method:`potdir_dict` invert keys from `paramfile_dict` to give 
    Returns: 
      {potential_file:potential_directory}
    """
    paramfile_dict = self.paramfile_dict()
    potdir = {}
    for k,v in paramfile_dict.items():
      potdir[v] = k
    return potdir

class Job(object):
  """
  class:'Job' collection of routines for generating 
  and submitting pbs submission scripts.
  """
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
    for target_dir in lst:
      target_dir = os.path.join(job_dir, target_dir)
      if regex == None:
        if os.path.isdir(target_dir) and target_dir != 'DFT':
          self.sub_pbs(dir, suffix=suffix, regex=regex)
        elif target_dir.split('_')[-1] == suffix:
          pbs_dir = os.path.join(sub_dir, target_dir)
          os.system("cd {0}; qsub fe{1}.pbs".format(pbs_dir, job_dir+'_'+suffix))
        else:
          pass
      else:
        if os.path.isdir(target_dir) and target_dir != 'DFT':
          self.sub_pbs(target_dir, suffix=suffix, regex=regex)
        elif regex.match(target_dir):
          try:
            target_dir  = '/'.join(target_dir.split('/')[:-1])
            name = target_dir.split('/')[-1]
            os.system("cd {0}; qsub fe{1}.pbs".format(target_dir, name))
          except:
            print 'Job Submit Failed'
        else:
          pass

class GBMaintenance(object):
  """
  :class'GBMaintenance' is a collection of maintenance routines for the GB database.
  Possible usages: regenerate all the csl lattices in the database
  or a subdirectory of the database, take a new grain boundary profile 
  picture for multiple directories, update the gb json information if a
  new grain boundary property is desired.
  """
  def __init__(self):
    self.materials = ['alphaFe']

  def retake_pic(self,fname, translate=False,toggle=False, confirm=True):
    """ 
    If using AtomEye take grain boundary profile pic in directory
    requires gb directory with gbid.xyz file in it. set confirm = False 
    to not prompt for overwrite.
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

  def remove_eo_files(self, path, num_deleted_files, dryrun=False):
    """
    Remove files with pattern matching jobname.[eo][0-9]+. 
    Returns: number of deleted files.
    """
    eo_regex = re.compile(r'[eo][0-9]+')
    lst = os.listdir(path)
    for filename in lst:
      filename = os.path.join(path, filename)
      if os.path.isdir(filename):
        rec_deleted_files = self.remove_eo_files(filename, 0)
        num_deleted_files += rec_deleted_files
      elif eo_regex.match(filename.split('.')[-1]):
        print filename
        os.remove(filename)
        num_deleted_files += 1
      else:
        pass
    return num_deleted_files

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
    new_json['H']    = cell[2,2]
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
	      self.fix_json(new_path)
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


class GBAnalysis(object):
  def __init__(self):
    pass

  def find_gb_json(self, path, j_list, filetype):
    """ 
    :method:`find_gb_json` Populates the list j_list with lists of the form
    [/directory_path/, /subgb_file_path].
    attributes:
      path     : root directory to begin recursive search
      j_list   : empty list  to populate
      filetype : 'subgb.json', 'gb.json'
    """
    lst = os.listdir(path)
    for filename in lst:
      filename = os.path.join(path, filename)
      if os.path.isdir(filename):
        self.find_gb_json(filename, j_list, filetype)
      elif filename.split('/')[-1] == filetype:
        j_list.append([path,filename])
      else:
        pass

  def extract_energies(self, material='alphaFe', or_axis='001', gb_type='tilt'):
    """
    :method:`extract_energies` pull GB formation energies using recursion.
    Go into a grain boundary directory, recurse through directory
    to find subgb.json files and gb_energies pull them out.
    Returns:
    list of dictionaries: each dictionary
    contains values gbid, orientation_axis, angle, boundary_plane, 
    param_file, and energies.
    """
    pot_param     = PotentialParameters()
    ener_per_atom = pot_param.gs_ener_per_atom()
    gb_files = []
#dir str lets us point to whatever material and orientation axisi it is we want in the database.
    dir_str  = os.path.join(material, or_axis)
    self.find_gb_json('{0}'.format(os.path.join(app.config['GRAIN_DATABASE'], dir_str)), gb_files, 'gb.json')
    grain_energies = []
    for grain in  gb_files:
      path = grain[0]
      with open(grain[1],'r') as f:
        j_dict    = json.load(f)
      subgb_files = []
      self.find_gb_json(path, subgb_files, 'subgb.json')
      calc_types  = []
#Find all calculation types associated with this grain
      for subgrain in subgb_files:
        with open(subgrain[1],'r') as f:
          try:
            sub_dict = json.load(f)
          except:
            pass
        try:
          if sub_dict['param_file'] not in calc_types:
            calc_types.append(sub_dict['param_file'])
        except KeyError:
          pass
#Initialize a dictionary of dictionaries for each calc type:
      gb_dict = {}
      for calc in calc_types:
        tmp_dict = {}
        tmp_dict['orientation_axis'] = j_dict['orientation_axis']
        tmp_dict['angle']            = j_dict['angle']*(180./np.pi)
        tmp_dict['boundary_plane']   = j_dict['boundary_plane']
        tmp_dict['energies']         = []
        tmp_dict['param_file']       = calc
        tmp_dict['gbid']             = j_dict['gbid']
        gb_dict[calc] = tmp_dict

      for subgrain in subgb_files:
        with open(subgrain[1],'r') as f:
          try:
            sub_dict = json.load(f)
            try:
              append_energy = True
              if 'iron_mish.xml' == sub_dict['param_file']:
                gb_ener = 16.02*((sub_dict['E_gb']-(ener_per_atom['iron_mish.xml']*float(sub_dict['n_at'])))/(2*sub_dict['A']))
              elif 'Fe_Mendelev.xml' == sub_dict['param_file']:  
                gb_ener = 16.02*((sub_dict['E_gb']-(ener_per_atom['Fe_Mendelev.xml']*float(sub_dict['n_at'])))/(2*sub_dict['A']))
              elif 'Fe_Ackland.xml' == sub_dict['param_file']:  
                gb_ener = 16.02*((sub_dict['E_gb']-(ener_per_atom['Fe_Ackland.xml']*float(sub_dict['n_at'])))/(2*sub_dict['A']))
              elif 'Fe_Dudarev.xml' == sub_dict['param_file']:  
                gb_ener = 16.02*((sub_dict['E_gb']-(ener_per_atom['Fe_Dudarev.xml']*float(sub_dict['n_at'])))/(2*sub_dict['A']))
              elif 'PotBH.xml' == sub_dict['param_file']:  
                gb_ener = 16.02*((sub_dict['E_gb']-(ener_per_atom['PotBH.xml']*float(sub_dict['n_at'])))/(2*sub_dict['A']))
              elif 'dft_vasp_pbe' == sub_dict['param_file']:  
                gb_ener = 16.02*((sub_dict['E_gb']-(ener_per_atom['dft_vasp_pbe']*float(sub_dict['n_at'])))/(2*sub_dict['A']))
              else:
                append_energy = False
                print 'Ground state energy not know for this potential!'
              if append_energy == True:
                gb_dict[sub_dict['param_file']]['energies'].append(gb_ener)
            except KeyError:
              pass
          except:
            pass
      for gdict in gb_dict.values():
        grain_energies.append(gdict)
    return grain_energies

  def calc_energy(self, gb_dict, param_file='PotBH.xml'):
    """
    :method:`calc_energy` given a subgb.json dictionary, and a potential
    calculate grainboundary energy.
    """
    pot_param     = PotentialParameters()
    ener_per_atom = pot_param.gs_ener_per_atom()
    try:
      gb_ener = 16.02*((gb_dict['E_gb']-(ener_per_atom[param_file]*float(gb_dict['n_at'])))/(2*gb_dict['A']))
    except KeyError:
      return None
    else:
      return gb_ener

  def pull_gamsurf(self, path="./",  potential="PotBH"):
    """
    :method:`pull_gamsurf` Loop over subgrain directories of a potential (default PotBH) 
    to find the minimum and maximum energies for the canonical grain, return
    a dictionary, with information about the lowest and highest energy structures.
    """
    potparams = PotentialParameters()
    paramfile_dict = potparams.paramfile_dict()
    subgb_files = []
    if os.path.isdir(os.path.join(path,potential)):
      self.find_gb_json(os.path.join(path,potential), subgb_files, 'subgb.json')
      gam_surfs   = []
      unconv      = []
#Only pulling for PotBH:
      for gb in subgb_files:
        with open(gb[1],'r') as f:
          gb_json = json.load(f)
        ener = self.calc_energy(gb_json, param_file=paramfile_dict[potential])
        if ener != None:
          gam_surfs.append((gb_json['rcut'], gb_json['rbt'][0], gb_json['rbt'][1], ener))
        else:
          unconv.append(gb[1])
      en_list     = [x[3] for x in gam_surfs]
      try:
        min_en      = min(en_list)
      except ValueError:
        return {'max_en':0.0, 'min_en':0.0, 'min_coords':[], 'max_coords':[]}
      else:
#Create lists of minimum energy structures (vx bxv rcut).
        min_coords  = [(gam[1], gam[2], gam[0]) for gam in filter(lambda x: round(x[3], 5) == round(min_en, 5), gam_surfs)]
        max_en      =  max(en_list)
        max_coords  = [(gam[1], gam[2], gam[0]) for gam in filter(lambda x: round(x[3], 5)==round(max_en, 5), gam_surfs)]
        gam_dict    = {'max_en':max_en, 'min_en':min_en, 'min_coords':min_coords, 'max_coords':max_coords}
    else:
      print "No potential directory:", potential, "found."
      gam_dict = {'max_en':0.0, 'min_en':0.0,'min_coords':[],'max_coords':[]}
    return gam_dict

  def plot_gamsurf(self, pot_dir='PotBH', rcut=None, print_gamsurf=False):
    """
    :method:`plot_gamsurf` return a dictionary for the gamma surface 
             and the directory with the lowest energy structure.
    """
    subgb_files = []
    self.find_gb_json(pot_dir, subgb_files, 'subgb.json')
    gam_surfs   = []
    for gb in subgb_files:
      with open(gb[1],'r') as f:
        gb_json = json.load(f)
      locen = self.calc_energy(gb_json)
      if locen is not None and locen >= 0.0:
        gam_surfs.append((gb_json['rcut'], gb_json['rbt'][0], gb_json['rbt'][1], locen, gb[0]))
      else:
        pass
    logging.info('gam_surfs', gam_surfs)
    if rcut is not None:
      gam_surf_rcut = filter(lambda x: x[0]==rcut, gam_surfs)
    else:
      gam_surf_rcut = gam_surfs

    if print_gamsurf:
      for count, gs in enumerate(gam_surf_rcut):
        print gs[1], gs[2], gs[3]
        if ((count+1)%6==0):
          print '\n'

    en_list = [x[3] for x in gam_surfs]
    en_list = [x for x in en_list if x is not None]
    min_en  = min(en_list)
    print 'Min Energy: ', min_en, 'J/m^{2}' 
    min_coords = filter(lambda x: round(x[3], 5) == round(min_en, 5), gam_surfs)
    print 'Coordinates of Min Energy Grain Boundaries:'
    for m in min_coords:
      print m
    max_en  = max(en_list)
    print 'Max Energy: ', max_en, 'J/m^{2}'
    print 'Coordinates of Max Energy Grain Boundaries:'
    max_coords = filter(lambda x: round(x[3], 5)==round(max_en, 5), gam_surfs)
    for m in max_coords:
      print m
    return min_coords, max_coords

  def list_all_unconverged(self, pattern='b111'):
    """
    :method:`list_all_unconverged` searches through subdirectories
    that  match against pattern for subgb.json files.
    Parameters:
      pattern: regex to pattern match against.
    Returns:
      Three lists converged_dirs, unconverged_dirs, dirs_missing_convergence_keys
    """
    jobdirs = glob.glob('{0}*'.format(pattern))
    jobdirs = filter(os.path.isdir, jobdirs)
    scratch = os.getcwd()
    converged_list   = []
    unconverged_list = []
    missing_key_list = []
    for job in jobdirs:
      os.chdir(os.path.join(scratch, job))
      subgb_files = []
      self.find_gb_json('./', subgb_files, 'subgb.json')
      for gb in subgb_files:
        with open(gb[1],'r') as f:
          gb_json = json.load(f)
        if 'converged' in gb_json:
          if gb_json['converged']:
            #print 'Converged', gb
            converged_list.append([job]+gb)
          elif not gb_json['converged']:
            #print 'Not Converged', gb
            unconverged_list.append([job]+gb)
        elif 'converged' not in gb:
            #print 'subgb missing converged key', gb
            missing_key_list.append(gb)
    os.chdir(scratch)
    return  converged_list, unconverged_list, missing_key_list

  def list_unconverged(self, pattern='001', potential='PotBH'):
    """
    :method:`list_unconverged` find all unconverged subgrains. 
    Convergence flag is in subgb.json file. 
    Parameters:
      pattern: Orientation axis
      potential: Interatomic Potential to search.
    Returns:
      Three lists converged_dirs, unconverged_dirs, dirs_missing_convergence_keys
    """
    jobdirs = glob.glob('{0}*'.format(pattern))
    jobdirs = filter(os.path.isdir, jobdirs)
    scratch = os.getcwd()
    converged_list   = []
    unconverged_list = []
    missing_key_list = []
    for job in jobdirs:
      os.chdir(os.path.join(scratch, job))
      subgb_files = []
      if os.path.isdir(potential):
        self.find_gb_json(potential, subgb_files, 'subgb.json')
        for gb in subgb_files:
          with open(gb[1],'r') as f:
            gb_json = json.load(f)
          if 'converged' in gb_json:
            if gb_json['converged']:
             # print 'Converged', gb
              converged_list.append([job]+gb)
            elif not gb_json['converged']:
             # print 'Not Converged', gb
              unconverged_list.append([job]+gb)
          elif 'converged' not in gb:
             # print 'subgb missing converged key', gb
              missing_key_list.append(gb)
    os.chdir(scratch)
    return  converged_list, unconverged_list, missing_key_list

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-e", "--extracten", action="store_true", help="Pull all energies for an orientation axis \
                                                    print lowest energies to terminal. List is ordered by angle.")
  parser.add_argument("-g", "--gam_min", action="store_true", help="Pull gamma surface for specified potential directory.")
  parser.add_argument("-d", "--directory", default="PotBH", help="Directory to search for min_en structure. (Default PotBH).")
  parser.add_argument("-m", "--material", help="The material we wish to query. Default (alphaFe).", default="alphaFe")
  parser.add_argument("-o", "--or_axis", help="Orientation axis.", default="001")
  parser.add_argument("-pt", "--potential", help="Potential file.", default ="PotBH.xml")
  args = parser.parse_args()
  analyze =  GBAnalysis()

  if args.extracten:
    or_axis = args.or_axis
    gb_list = analyze.extract_energies(or_axis=or_axis)
    for gb in sorted(gb_list, key = lambda x: x['angle']):
      if gb['param_file'] == args.potential:
        try:
          print gb['param_file'], gb['angle'], min(gb['energies']), max(gb['energies'])
        except ValueError:
          print 'No Valid Energy: ', gb['param_file'], gb['angle']
  
  if args.gam_min:
#   Search potential directory for the gamma surface it contains
#   for all the cutoff radii.
    analyze.gam_min(directory=args.directory)
