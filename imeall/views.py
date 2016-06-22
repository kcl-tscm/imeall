# -*- coding: utf-8 -*-
import os
import re
import json
from flask   import Flask, request, session, g, redirect, url_for, abort,\
                    render_template, flash, send_file, jsonify, make_response
from imeall  import app
from imeall.models import GBAnalysis
# Unique key is BBBAAAACCC
# Common axis=[BBB], misorientation angle=AAAA, and GB plane = (CCC).
# temporary table should be replaced by database.
# Associated with each Grain Boundary we need a (many) calculation(s) object(s)
# these should have some tag i.e. DFT-VASP-PBE. So clearly a one to many
# relationship.
# One GrainBoundary table with grain_boundary id (unique tag), name, raw atom structure
# coincident site lattice etc, and a bunch of pointers to a CalculationTable
# Each CalculationTable will have its own unique tag (DFT-VASP-PBE) along with the properties
# Calculate: total energy, forces, atoms, magnetic moments.
# Table energies should be populated in eV:
# Each calculation should have the atoms object attached to it.
grain_boundaries = {}
calculations      = {}
grain_boundaries['0000000000'] = {'title': 'Ideal Crystal ',  'gb_id':'0000000000'}
grain_boundaries['1107053111'] = {'title': 'Sigma(3)  (111)', 'gb_id':'1107053111'}
grain_boundaries['1105048332'] = {'title': 'Sigma(11) (332)', 'gb_id':'1105048332'}
grain_boundaries['1106000112'] = {'title': 'Sigma(3) (112)',  'gb_id':'1106000112'}
calculations['0000000000'] = {'VASP-DFT-PBE' : {'E0':-8.23807, 'DFT-mag': 2.2238, 'nat':1}, 'IP-EAM-MISH':{'E0': -4.2701, 'nat':1}}
calculations['1107053111'] = {'VASP-DFT-PBE' : {'E0':-406.154623782, 'nat':96, 'A': 27.7436434255}}
calculations['1105048332'] = {'IP-EAM-MISH'  : {'E0':-382.847802363, 'nat':90, 'A':18.7825353894 }}
calculations['1106000112'] = {'IP-EAM-MISH'  : {'E0':-196.171, 'nat':46, 'A': 9.80885920049}}
#files to display in browser:
valid_extensions = ['xyz', 'json']
vasp_files       = ['IBZKPT', 'INCAR', 'CHG', 'CHGCAR', 'DOSCAR', 'EIGENVAL', 
                    'KPOINTS', 'OSZICAR', 'OUTCAR', 'PCDAT', 'POSCAR',
                    'POTCAR', 'WAVECAR', 'XDATCAR']
# Currently the database connection is just a path name to our grain boundary
# database stored as a file tree. I don't necessarily see any reason not to exploit
# the existing filesystem and tools associated for searching. Why?
#   1) The nature of the work pattern is I'll want to be able to rummage around
#      in the different grain directories and subgrain directories to run quippy
#      scripts etc, or copy subgraindirs with DFT stuff in them.
#      Any work generated during this should just reside in the
#      directory and I will only make what I want visible to the imeall browser
#      xyz files for structures/forces, POSCAR, svg files for images, and json
#      the rest will be invisible.
#   2) By storing everything in the file system and only making what I want
#      visible we still preserve a hierarchical relationship for all the grain
#      boundaries within a class of materials and for different classes of
#      materials. This tree structure would resemble a collection of documents and
#      could be mapped onto a NoSQL type of database fairly easily. 
@app.before_request
def before_request():
  g.gb_dir = app.config['GRAIN_DATABASE']

@app.route('/')
def home_page():
  materials = os.listdir(g.gb_dir)
  return render_template('imeall.html', materials=materials)

@app.route('/<material>/')
def material(material):
  path     = os.path.join(g.gb_dir, material)
  url_path = material
  orientations = []
  files =  os.listdir(path)
  for filename in files: 
    tmp_file = os.path.join(path, filename)
    if os.path.isdir(tmp_file):
      orientations.append(filename)
  return render_template('material.html', url_path=url_path, orientations=orientations)

@app.route('/analysis/')
def analysis():
  '''
    This view collates data from the grainboundary database.
  '''
# User chooses what orientation angle to look at via a GET argument:
  or_axis = request.args.get('or_axis', '001')
  analyze = GBAnalysis()
  gb_list = analyze.extract_energies(or_axis=or_axis)
  gbdat   = []
# Creates list of grain boundaries ordered by angle.
  for gb in sorted(gb_list, key = lambda x: x['angle']):
    try:
      min_en = min([x for x in gb['energies'] if x > 0.])
      max_en = max(gb['energies'])
      try:
        gbdat.append({'param_file': gb['param_file'], 'or_axis':' '.join(map(str, gb['orientation_axis'])), 
                      'angle':gb['angle'], 'min_en':min_en, 
                      'max_en':max_en,
                      'bp':' '.join(map(str, map(int, gb['boundary_plane'])))})
      except KeyError:
        pass
    except ValueError:
      pass
  gbdat = sorted(gbdat, key=lambda x: x['angle'])
  return render_template('analysis.html', gbdat=json.dumps(gbdat))

@app.route('/orientation/<path:url_path>/<orientation>/')
def orientations(url_path, orientation):
  url_path = url_path+'/'+orientation
  path     = os.path.join(g.gb_dir, url_path)
  grains    = []
#Only valid directories beginning with orientation axis will be shown.
  for thing in os.listdir(path):
    if thing[:3] == orientation: 
      grains.append(thing.strip()) 
  return render_template('orientation.html', url_path=url_path, grains=grains)

def make_tree(path):
  tree = dict(name=os.path.basename(path), children=[], fullpath='')
  try: 
    lst = os.listdir(path)
    lst = map(os.path.basename, lst)
  except OSError:
    pass  
  else:
    for name in lst:
      filename = os.path.join(path,name)
      if os.path.isdir(filename):
        tree['children'].append(make_tree(filename))
      else:
#append file if it is a relevantwith its route:
        extension = name.split(".")[-1]
        if name in vasp_files or extension in valid_extensions:
          tree['children'].append(dict(name=name, fullpath=filename))
  return tree

def extract_json(path, json_files):
  lst = os.listdir(path)
  for f in lst:
    f = os.path.join(path,f)
    if os.path.isdir(f):
      extract_json(f, json_files)
    else:
      if f.split(".")[-1] == 'json' and (f.split("/")[-1])[:3]=='sub':
        json_files.append(f)
      else:
        pass

@app.route('/grain/<path:url_path>/<gbid>/')
def grain_boundary(url_path, gbid):
  url_path  = url_path+'/'+gbid
  path = os.path.join(g.gb_dir, url_path)
  with open(os.path.join(path, 'gb.json'),'r') as json_file:
    gb_info = json.load(json_file)
  stuff = []
  tree = make_tree(path)
# Get all the subdirectory json files so that
# we can easily compare normalized grain boundary formation
# energies, hydrogen interstitials etc.
  json_files = []
  extract_json(path, json_files)
  subgrains = []
  for i, path in enumerate(json_files):
    try: 
      subgrains.append([json.load(open(path,'r')), i])
    except:
      pass
  print 'PATH', path, 'stuff'
  print 'URL_PATH', url_path 
  return render_template('grain_boundary.html', gbid=gbid, url_path=url_path, 
                          stuff=stuff, gb_info=gb_info, tree=tree,
                          subgrains=subgrains)

#Check for Ovito in different paths.
def run_ovito(target_dir, filename):
  """ run_ovito is meant to launch the ovito viewer application with the
      associated grainboundary trajectory file loaded, the os command should
      ensure we are in the working directory so that any modifications, or
      videos generated will be saved in the correct place.
  """
  ovito = "~/ovito-2.6.1-x86_64/bin/ovito"
  target_dir = '/'.join(target_dir.split('/')[:-1])
  print target_dir
  print filename
  if os.path.isfile(ovito):
    os.system("cd {0}; {1} {2}&".format(target_dir, ovito, filename))
    variable = raw_input('Continue?')
  elif os.path.isfile('/Users/lambert/Ovito.app/Contents/MacOS/ovito'):
    os.system('cd {0}; /Users/lambert/Ovito.app/Contents/MacOS/ovito {1}&'.format(target_dir, filename))
  else: 
    return flash('Ovito not installed in one of the usual locations.')

#This route serves images from the grain boundary directory.
@app.route('/img/<path:filename>/<gbid>/<img_type>')
def serve_img(filename, gbid, img_type):
#should consider security stuff here as well... flask.safe_join()
  print 'FILENAME', filename
  print 'GBID', gbid
  img = os.path.join(filename,'{0}.png'.format(gbid))
  if img_type =='struct':
    img  = app.config['GRAIN_DATABASE']+'/'+filename+'/{0}.png'.format(gbid)
  elif img_type =='csl':
    img  = app.config['GRAIN_DATABASE']+'/'+filename+'/csl_{0}.svg'.format(gbid)
  elif img_type =='pot':
    pot_dir = '/Users/lambert/pymodules/imeall/imeall/potentials'
    img     = pot_dir+'/'+filename
  return send_file(img)

@app.route('/textfile/<gbid>/<path:filename>')
def serve_file(gbid, filename):
#Serves text files should we want to inspect these in the browser:
  textpath = request.args.get('textpath')
  with open('{0}'.format(textpath),'r') as text_file:
    text = text_file.read()
  if filename.split(".")[-1] == 'xyz':
    text = text.split('\n')
    run_ovito(textpath, filename)
    return render_template('text.html', text=text)
  elif filename == 'OSZICAR':  
    rmm_regex  = re.compile(r'RMM:\s+([-+0-9.E\s]+)')
    osz_list = rmm_regex.findall(text)
    osz = []
    line     = {}
    for i in range(len(osz_list[:80])):
      split_list = map(float, osz_list[i].split())
      if i + 1  == split_list[0]:
        osz.append({'N'     : split_list [0],
                    'E'     : split_list [1],
                    'dE'    : split_list [2],
                    'd_eps' : split_list [3],
                    'ncg'   : split_list [4],
                    'rms'   : split_list [5]})
      else:
        break
    return render_template('oszicar.html', osz=json.dumps(osz))
  else:
    text = text.split('\n')
    return render_template('text.html', text=text)

@app.route('/eam_pot/<path:filename>')
def eam_pot(filename):
#found this gist at https://gist.github.com/wilsaj/862153
  import random
  import datetime
  import StringIO
  import xml.etree.ElementTree as ET
  import matplotlib.pyplot as plt

  from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
  from matplotlib.figure import Figure
  from matplotlib.dates import DateFormatter
  

  fig = Figure()
  ax = fig.add_subplot(111)
  FVR = []
  pot_dir  = '/Users/lambert/pymodules/imeall/imeall/potentials'
  pot_xml  = pot_dir+'/'+filename
  with open(pot_xml, 'r') as f:
    tree = ET.parse(f)
  root = tree.getroot()
  R = [float(x.attrib['y']) for x in root[0][0]]
  F = [float(x.attrib['y']) for x in root[0][1]]
  V = [float(x.attrib['y']) for x in root[1][0]]
  xpts = []
  for x in root.findall('./per_type_data/spline_rho/'):
    xpts.append(x.attrib['r'])
  rho = [float(x.attrib['r']) for x in root[0][1]]
  fig, ax = plt.subplots(3,1)
  ax[0].set_ylim([-0.5,5])
  ax[0].set_xlim([1.,5.5])
  ax[1].set_xlim([1.,5.5])
  ax[0].set_ylabel('Density')
  ax[0].plot(xpts, R)
  ax[1].set_ylim([-5,50])
  ax[1].set_ylabel('V')
  ax[1].set_xlabel('Radius (A)')
  ax[1].plot(xpts,V)
  ax[2].set_ylim([-50,20])
  ax[2].set_ylabel('Embedding Term')
  ax[2].set_xlabel('Density (arb. units)')
  ax[2].plot(rho, F)
  canvas     = FigureCanvas(fig)
  png_output = StringIO.StringIO()
  canvas.print_png(png_output)
  response   = make_response(png_output.getvalue())
  response.headers['Content-Type']= 'image/png'
  return response

