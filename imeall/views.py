# -*- coding: utf-8 -*-
import os
import re
from imeall  import app
from flask   import Flask, request, session, g, redirect, url_for, abort,\
                    render_template, flash, send_file, jsonify
import json
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


valid_extensions = ['xyz']
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
  orientations     = os.listdir(path)
  return render_template('material.html', url_path=url_path, orientations=orientations)

@app.route('/orientation/<path:url_path>/<orientation>/')
def orientations(url_path, orientation):
  url_path = url_path+'/'+orientation
  path     = os.path.join(g.gb_dir, url_path)
  print 'Orientation path', path
  print  url_path
  grains    = []
#Only valid directories beginning with orientation axis will be shown.
  for thing in os.listdir(path):
    if thing[:3] == orientation: 
      grains.append(thing.strip()) 
      print thing.strip()
  return render_template('orientation.html', url_path=url_path, grains=grains)

#This seems like a very useful pattern to learn!
def make_tree(path):
  tree = dict(name=os.path.basename(path), children=[], fullpath='')
  try: 
    lst = os.listdir(path)
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

@app.route('/grain/<path:url_path>/<gbid>/')
def grain_boundary(url_path, gbid):
  url_path  = url_path+'/'+gbid
  path = os.path.join(g.gb_dir, url_path)
  with open(os.path.join(path, 'gb.json'),'r') as json_file:
    gb_info = json.load(json_file)
  stuff = []
# really  what we want to do here is to 
# walk through the subdirectories
#  for thing in os.listdir(path):
#    if os.path.isdir(os.path.join(path,thing)):
#      stuff.append(thing) 
  tree = make_tree(path)
  print 'PATH', path, 'stuff'
  print 'URL_PATH', url_path 
  return render_template('grain_boundary.html', gbid=gbid, url_path=url_path, 
                          stuff=stuff, gb_info=gb_info, tree=tree)

#Check for Ovito in different paths.
@app.route('/ovito/<path:target_dir>/<gbid>/<input_type>')
def run_ovito(target_dir, gbid, input_type):
  """ run_ovito is meant to launch the ovito viewer application with the
      associated grainboundary trajectory file loaded, the os command should
      ensure we are in the working directory so that any modifications, or
      videos generated will be saved in the correct place.
  """
  ovito = "~/ovito-2.6.1-x86_64/bin/ovito"
  if os.path.isfile(ovito):
    os.system("cd {0}; {1} {2}.xyz".format(target_dir, ovito, os.path.join(target_dir, gbid)))
    variable = raw_input('Continue?')
  elif os.path.isfile('/Users/lambert/Ovito.app/Contents/MacOS/ovito'):
    os.system('cd {0}; /Users/lambert/Ovito.app/Contents/MacOS/ovito {0}.xyz'.format(os.path.join(target_dir, gbid)))
  else: 
    return flash('Ovito not installed in one of the usual locations.')

#This route serves images from the grain boundary directory.
@app.route('/img/<path:filename>/<gbid>/<img_type>')
def serve_img(filename, gbid, img_type):
#should consider security stuff here as well... flas.safe_join()
  print 'FILENAME', filename
  print 'GBID', gbid
  img  = os.path.join(filename,'{0}.png'.format(gbid))
  if img_type =='struct':
    img  = app.config['GRAIN_DATABASE']+'/'+filename+'/{0}.png'.format(gbid)
  elif img_type =='csl':
    img  = app.config['GRAIN_DATABASE']+'/'+filename+'/csl_{0}.svg'.format(gbid)
#Might want to use flask.send_from_directory() here.
  return send_file(img)

@app.route('/textfile/<path:textpath>')
def serve_file(textpath):
#serves txt files should we want to inspect these in the browser.
  with open('{0}'.format(textpath),'r') as text_file:
    text = text_file.read()
  return 'Hello World'


