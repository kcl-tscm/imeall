# -*- coding: utf-8 -*-
import os
import re
import json
import logging
import subprocess
from imeall import app
from flask import Flask, request, session, g, redirect, url_for, abort,\
                  render_template, flash, send_file, jsonify, make_response, safe_join
from imeall.models import GBAnalysis
from gb_models import serialize_vector, GRAIN_DATABASE, DATABASE, GrainBoundary, SubGrainBoundary,\
                      deserialize_vector_int
from models import PotentialParameters

calculations      = {}
calculations['0000000000'] = {'VASP-DFT-PBE' : {'E0':-8.23807, 'DFT-mag': 2.2238, 'nat':1}, 'IP-EAM-MISH':{'E0': -4.2701, 'nat':1}}
calculations['1107053111'] = {'VASP-DFT-PBE' : {'E0':-406.154623782, 'nat':96, 'A': 27.7436434255}}
calculations['1105048332'] = {'IP-EAM-MISH'  : {'E0':-382.847802363, 'nat':90, 'A':18.7825353894 }}
calculations['1106000112'] = {'IP-EAM-MISH'  : {'E0':-196.171, 'nat':46, 'A': 9.80885920049}}

#files that Imeall server can wants to display in browser:
valid_extensions = ['xyz', 'json', 'mp4', 'png','day']
vasp_files       = ['IBZKPT', 'INCAR', 'CHG', 'CHGCAR', 'DOSCAR', 'EIGENVAL', 
                    'KPOINTS', 'OSZICAR', 'OUTCAR', 'PCDAT', 'POSCAR',
                    'POTCAR', 'WAVECAR', 'XDATCAR']

@app.before_request
def before_request():
  g.gb_dir = app.config['GRAIN_DATABASE']

@app.route('/')
def home_page():
  """Overview of imeall database. Links to material views.
  """
  materials = app.config["MATERIALS"]
  return render_template('imeall.html', materials=materials)

@app.route('/<material>/')
def material(material):
  """View of available orientation axes for a particular material.
  """
  path         = os.path.join(app.config['GRAIN_DATABASE'], material)
  url_path     = material
  orientations = []
  files =  os.listdir(path)
  files.sort()
  for filename in files: 
    tmp_file = os.path.join(path, filename)
    if os.path.isdir(tmp_file):
      orientations.append(filename)
  return render_template('material.html', url_path=url_path, orientations=orientations)

@app.route('/orientation/<path:url_path>/<orientation>/')
def orientations(url_path, orientation):
  """View to list different orientation axes in the material database.
  """

#can only handle three digit or_axis atm.
  url_path = url_path+'/'+orientation
  path     = os.path.join(g.gb_dir, url_path)

#load serialized grain data
  with open(os.path.join(path, 'or_axis.json'), 'r') as json_file:
    oraxis = json.load(json_file)
  oraxis = oraxis['oraxis']
  gb_type = request.args.get('gb_type', 'tilt')

  if gb_type == 'tilt':
    gbs   = (GrainBoundary.select().where(GrainBoundary.orientation_axis==oraxis)
                          .where(GrainBoundary.boundary_plane !=oraxis)
                          .order_by(GrainBoundary.angle))
  elif gb_type == 'twist':
    gbs   = (GrainBoundary.select().where(GrainBoundary.orientation_axis==oraxis)
                          .where(GrainBoundary.boundary_plane == oraxis)
                          .order_by(GrainBoundary.angle))
  else:
    gbs   = GrainBoundary.select().where(GrainBoundary.orientation_axis==oraxis).order_by(GrainBoundary.angle)

#only valid directories beginning with orientation axis will be shown.
  grains = []
  for gb in gbs:
    grains.append({'gbid':gb.gbid, 'angle':round(gb.angle*180/(3.14159), 2), 'bp':deserialize_vector_int(gb.boundary_plane)}) 
  return render_template('orientation.html', url_path=url_path, grains=grains)

@app.route('/grain/<path:url_path>/<gbid>/')
def grain_boundary(url_path, gbid):
  """Top view for a canonical grain boundary. CSL 
  lattice, and list of subgrain directories, and energies.

  Todo: 
    Replace make_tree with an SQL query.
  """

  url_path  = url_path+'/'+gbid
  path      = os.path.join(app.config["GRAIN_DATABASE"], url_path)
  with open(os.path.join(path, 'gb.json'),'r') as json_file:
    gb_info = json.load(json_file)
  stuff = []
  tree  = make_tree(path)
  json_files = []
  extract_json(path, json_files)
  subgrains  = []  
  subgrainsj = []
  for i, gb_path in enumerate(json_files):
    subgrains.append([json.load(open(gb_path, 'r')), i])
    subgrainsj.append(json.load(open(gb_path, 'r')))

#find lowest energy structure for every potential in the database
  analyze  = GBAnalysis()
  potparams = PotentialParameters()
  paramfile_dict = potparams.paramfile_dict()
  gam_dict = {}
  for potdir in paramfile_dict.keys():
    gam_dict[potdir] = analyze.pull_gamsurf(path=path, potential=potdir) 
  return render_template('grain_boundary.html', gbid=gbid, url_path=url_path,
                          gb_info=gb_info, flare_root=json.dumps(tree), subgrains=subgrains, 
                          subgrainsj=json.dumps(subgrainsj), gam_dict=gam_dict)

@app.route('/analysis/')
def analysis():
  """This view collates data from the grain boundary database
  and forwards it to d3 analysis tools.
  """

# User chooses what orientation angle to look at via a GET argument:
# This should be a separate Table.
  pot_param = PotentialParameters()
  ener_per_atom = pot_param.gs_ener_per_atom()
  or_axis = request.args.get('oraxisselect', default='001')
  gb_type = request.args.get('gbtypeselect', default='tilt')
  material = request.args.get('materialselect', default='alphaFe')
  gbdat = []
  oraxis = ','.join([c for c in or_axis])
# Creates list of grain boundaries ordered by angle.
  for potential in ener_per_atom.keys():
# GrainBoundary Energies in J/m^{2}
    if gb_type == 'tilt':
      gbs = (GrainBoundary.select()
                         .where(GrainBoundary.orientation_axis==oraxis)
                         .where(GrainBoundary.boundary_plane != oraxis))
    elif gb_type == 'twist':
      gbs = (GrainBoundary.select()
                         .where(GrainBoundary.orientation_axis==oraxis)
                         .where(GrainBoundary.boundary_plane == oraxis))
    else:
      sys.exit('Invalid gb_type!')
    for gb in gbs.order_by(GrainBoundary.angle):
      subgbs = (gb.subgrains.select(GrainBoundary, SubGrainBoundary)
                      .where(SubGrainBoundary.potential==potential)
                      .join(GrainBoundary)
                      .order_by(SubGrainBoundary.E_gb)
                      .dicts())
      logging.debug(gb.gbid)
      subgbs = [(16.02*(subgb['E_gb']-float(subgb['n_at']*ener_per_atom[potential]))/(2.0*subgb['area']), subgb) for subgb in subgbs]
      subgbs.sort(key = lambda x: x[0])
      if (len(subgbs) > 0):
        if (subgbs[0][1]['converged'] == True and subgbs[0][0] < 3.0):
          gbdat.append({'param_file' : potential,
                      'or_axis'    : ' '.join(map(str, subgbs[0][1]['orientation_axis'].split(','))),
                      'angle'      : subgbs[0][1]['angle']*(180./(3.14159)),
                      'min_en'     : subgbs[0][0],
                      'bp'         : ' '.join(map(str, map(int, deserialize_vector_int(subgbs[0][1]['boundary_plane'])))),
                      'url'        : 'http://127.0.0.1:5000/grain/alphaFe/'
                                    + ''.join(map(str, deserialize_vector_int(subgbs[0][1]['orientation_axis'])))
                                    + '/' + gb.gbid})
      else:
        pass
  return render_template('analysis.html', gbdat=json.dumps(gbdat))

def make_tree(path):
  """Recurse through subgrain directories collecting json and png files.
  """

  tree = dict(name=os.path.basename(path), children=[], fullpath='')
  try: 
    lst = os.listdir(path)
    lst = map(os.path.basename, lst)
  except OSError:
    pass  
  else:
    for name in lst:
      filename = os.path.join(path, name)
      if os.path.isdir(filename):
        tree['children'].append(make_tree(filename))
      else:
#append file if it is a relevant with its route:
        extension = name.split(".")[-1]
        #print filename
        filename = os.path.relpath(filename, app.root_path)
        if (name in vasp_files) or (extension in valid_extensions):
          tree['children'].append(dict(name=name, fullpath = url_for('serve_struct', filename='apathtoafile', textpath=filename)))
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

#check for Ovito in different paths.
def run_ovito(filename):
  """Launches the Ovito application with the
  associated grain boundary trajectory file loaded. 
  The alias 'ovito' must be set in environment. 
  This function can only be used if running server 
  locally with a local copy of the database.
  """

  try:
    ovito = os.environ["OVITO"]
  except KeyError:
    flash('No path to OVITO found in the environment')
  #app.logger.debug(ovito, os.path.dirname(filename),  os.path.basename(filename))
  job = subprocess.Popen("{0} {1}".format(ovito, os.path.basename(filename)).split(), cwd=os.path.dirname(filename))


@app.route('/struct/<path:filename>/<path:textpath>')
def serve_struct(filename, textpath=None):
  """View for serving structure files to clients.
  """

  if textpath.endswith('xyz'):
    if app.config["RUN_OVITO"]:
      run_ovito(os.path.join(app.config['GRAIN_DATABASE'], textpath))
      flash('running ovito')
      return redirect(request.referrer)
    else:
      return send_file(textpath)

  elif textpath.endswith('json'):
    with open(os.path.join(app.config['GRAIN_DATABASE'], textpath),'r') as f:
      json_dict = json.load(f)
    return render_template('json_dict.html', json_dict=json_dict)
  elif textpath.endswith('png'):
    print os.path.join(app.config['GRAIN_DATABASE'], textpath)
    return send_file(os.path.join(app.config['GRAIN_DATABASE'], textpath))
  else:
    return redirect(request.referrer)
  

@app.route('/img/<path:filename>/<gbid>/<img_type>')
def serve_img(filename, gbid, img_type):
  """Serve image_file to the browser.
  """

  img = os.path.join(filename,'{0}.png'.format(gbid))
  if img_type =='struct':
    img = app.config['GRAIN_DATABASE']+'/'+filename+'/{0}.png'.format(gbid)
  elif img_type =='csl':
    img = app.config['GRAIN_DATABASE']+'/'+filename+'/csl_{0}.svg'.format(gbid)
  elif img_type =='pot':
    pot_dir = os.path.join(app.root_path, 'potentials')
    img = pot_dir+'/'+filename
  elif img_type == 'gen':
    img = app.config['GRAIN_DATABASE']+'/'+filename
  else:
    img = 'NO IMAGE'
  return send_file(img)

@app.route('/servefile/')
def serve_file(textpath):
  """Serve different common file types to the browser.
  """

  with open('{0}'.format(textpath), 'r') as text_file:
    text = text_file.read()
  if textpath.endswith('xyz'):
    if app.config["RUN_OVITO"]:
      run_ovito(textpath)
      return render_template('text.html', text=text)
    else:
      return send_file(textpath)
  elif textpath.endswith('json'):
    with open(textpath, 'r') as f:
      j_file = json.load(f)
    return render_template("render_json.html", j_file=j_file)
  elif textpath.endswith('png'):
    return send_file(textpath)
  elif textpath == 'OSZICAR':  
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
  elif textpath.endswith('mp4'):
    subprocess.Popen('open {0}'.format(textpath).split()) 
    return render_template('text.html', text='Playing Video')
  else:
    text = text.split('\n')
    return render_template('text.html', text=text)

@app.route('/eam_pot/<path:filename>')
def eam_pot(filename):
  """Uses matplotlib to inspect xml potential files in the database.
  """

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
  pot_dir  = os.path.join(app.root, 'potentials')
  pot_xml  = pot_dir+'/'+filename
#Based on gist at https://gist.github.com/wilsaj/862153.
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

