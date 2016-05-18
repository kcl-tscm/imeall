import os
import sys
from flask  import Flask, request, session, g, redirect
from  flask import    url_for, abort, render_template, flash
from   quippy import Atoms
from   imeall import app
import imeall.slabmaker.slabmaker as slabmaker
import json
from imeall.models import GBAnalysis

class TestDB (object):
  '''
  Logic and consistency tests for the Grainboundary database.
  in particular routines to test all grain and subgrain json files
  are well formed (i.e.) contain logical keys, sigma indices are all
  integers, standard naming patterns for suffices.
  '''
  def __init(self):
    pass

  def extract_json(self, path, json_files):
    lst = os.listdir(path)
    for f in lst:
      f = os.path.join(path,f)
      if os.path.isdir(f):
        self.extract_json(f, json_files)
      else:
        if f.split(".")[-1] == 'json':
          json_files.append(f)
        else:
          pass


if __name__ == '__main__':
  j_files = []
  db_test = TestDB()
  gb_extract = GBAnalysis()
  j_files    = []
  gb_extract.find_gb_json('./110', j_files,'subgb.json')
  for j_file in j_files:
    j_dict = json.load(open(j_file[1],'r'))
    if 'n_at' not in j_dict.keys():
      print j_file, 'POORLY FORMED'
      with open(j_file[1] ,'r') as f:
        json_file = json.load(f)
      print j_file
      gbid              = j_file[0].split('/')[-1]
      print j_file, gbid
      try:
        at = Atoms('{0}.xyz'.format(os.path.join(j_file[0],gbid.split('_')[0]+'_n12d2.0')))
      except:
        try:
          at = Atoms('{0}.xyz'.format(os.path.join(j_file[0],gbid.split('_')[0]+'_n0d2.0')))
        except:
          try:
            at = Atoms('{0}.xyz'.format(os.path.join(j_file[0],gbid.split('_')[0]+'_n16d2.0')))
          except:
            try:
              at = Atoms('{0}.xyz'.format(os.path.join(j_file[0],gbid.split('_')[0]+'_n16d2.0')))
            except:
              try:
                at = Atoms('{0}.xyz'.format(os.path.join(j_file[0],gbid.split('_')[0]+'_n8d2.0')))
              except:
                print 'FAILED'
                break
      json_file['n_at'] = len(at)
      json_file['gbid'] = gbid
      with open(j_file[1] ,'w') as f:
        json.dump(json_file, f)
    if 'n_unit_cell' in j_dict.keys():
      print j_file, 'HAS n_unit_cell'
      db_test.update_json(j_file)
