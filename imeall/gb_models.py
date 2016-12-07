import os
import sys
import argparse
import json
from   peewee import *
from   datetime import datetime, timedelta
from   models import GBAnalysis

GRAIN_DATABASE = "/home/lambert/pymodules/imeall/imeall/grain_boundaries/"
DATABASE       = "gb_database.sql"
database       = SqliteDatabase(DATABASE)
class BaseModel(Model):
  class Meta():
    database = database

class GrainBoundary(BaseModel):
  """
  Canonical Parent Grain
  Vectors are serialized to csv. We may want to 
  separate this out into tables if heavy searching
  becomes necessary.
  :params: angle misorientation angle in radians.
  """
  gb_type          = CharField()
  boundary_plane   = CharField()
  z_planes         = CharField()
  orientation_axis = CharField()
  n_at             = IntegerField()
  coincident_sites = IntegerField()
  angle            = FloatField()
  height           = FloatField()
  area             = FloatField()
  notes            = TextField(default="")
  path             = CharField()
  gbid             = CharField(unique=True)

class SubGrainBoundary(BaseModel):
  """
  :path: relative to the grainboundary database root.
  :params: rbt rigid body translations.
  :params: grain_boundary every grain is a subgrain of the GrainBoundary Class.
  """
  grain_boundary = ForeignKeyField(GrainBoundary, "subgrains")
  converged      = BooleanField()
  rbt            = CharField()
  path           = CharField()
  potential      = CharField()
  rcut           = FloatField()
  area           = FloatField()
  n_at           = IntegerField()
  E_gb           = FloatField(default=0.0)
  E_gb_init      = FloatField(default=0.0)
  notes          = TextField(default="")
  gbid           = CharField(unique=True)

class Fracture(BaseModel):
  """
  :params:G stress energy release rate.
  :params:strain_rate.
  :params:sim_T simulation temperature.
  """
  fracture_system = CharField()
  G               = FloatField()
  strain_rate     = FloatField()
  sim_T           = FloatField()
  notes           = CharField(default="")

class Dislocation(BaseModel):
  gbid  =  CharField()

def serialize_vector(vector):
  return ','.join(map(str, vector))

def deserialize_vector_float(ser_vec):
  return map(float, servec.split(','))

def deserialize_vector_int(ser_vec):
  return map(int, servec.split(','))

def create_tables(database):
  """
  :method:`create_tables` 
  """
  database.create_tables([GrainBoundary,SubGrainBoundary])

def populate_db(or_axis='001'):
  analyze  = GBAnalysis()
  dir_str  = os.path.join('alphaFe', or_axis)
  gb_files = []
#grab list of gb.json files.
  analyze.find_gb_json('{0}'.format(os.path.join(GRAIN_DATABASE, dir_str)), gb_files, 'gb.json')
  for gb in gb_files:
    print gb[0], gb[1]
    with open(gb[1], 'r') as f:
      gb_json = json.load(f)
    print gb_json
    gb_dict = {"gb_type"          : gb_json['type'],
               "n_at"             : gb_json['n_at'],
               "boundary_plane"   : serialize_vector(map(int, gb_json['boundary_plane'])),
               "orientation_axis" : serialize_vector(map(int, gb_json['orientation_axis'])),
               "z_planes"         : serialize_vector(gb_json['zplanes']),
               "coincident_sites" : gb_json['coincident_sites'],
               "angle"            : gb_json['angle'],
               "height"           : gb_json['H'],
               "area"             : gb_json['A'],
               "notes"            : "",
               "path"             : os.path.relpath(gb[0], "/home/lambert/pymodules/imeall/imeall/grain_boundaries/"),
               "gbid"             : gb_json['gbid']
              }
#   GB_model_object = GrainBoundary.create(**gb_dict)
    GB_model_object = {}
    subgb_files = []
    analyze.find_gb_json('{0}'.format(gb[0]), subgb_files, 'subgb.json')
    for subgb in subgb_files:
      print 'SUBGB', subgb
      with open(subgb[1],'r') as f:
        subgb_json = json.load(f)
      try: 
        converged = subgb_json['converged']
      except KeyError:
        converged = False
      try:
        E_gb = subgb_json["E_gb"]
      except KeyError:
        E_gb = 0.0
      try:
        E_gb_init=subgb_json["E_gb_init"]
      except KeyError:
        E_gb_init = 0.0
      try:
        gbid = subgb_json["gbid"]
      except KeyError:
        gbid = subgb_json["name"]

      subgb_dict = {"grain_boundary"   : gb_model_object,
                    "converged"        : converged,
                    "E_gb_init"        : E_gb_init, 
                    "potential"        : subgb_json["param_file"],
                    "rbt"              : serialize_vector(subgb_json['rbt']),
                    "path"             : os.path.relpath(subgb[0], "/home/lambert/pymodules/imeall/imeall/grain_boundaries/"),
                    "area"             : gb_json['A'],
                    "rcut"             : subgb_json["rcut"],
                    "n_at"             : gb_json['n_at'],
                    "E_gb"             : E_gb,
                    "notes"            : "",
                    "gbid"             : gbid
                  }
#      SubGrainBoundary.create(**subgb_dict)        
if __name__=="__main__":
  populate_db("110")

