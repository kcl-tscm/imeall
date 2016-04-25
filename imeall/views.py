# -*- coding: utf-8 -*-
import os
import re
import peewee
from imeall  import app
from flask   import Flask, request, session, g, redirect, url_for, abort, render_template, flash


#
# Unique key is BBBAAAACCC
# Common axis=[BBB], misorientation angle=AAAA, and GB plane = (CCC).
# temporary table should be replaced by database.
# Associated with each Grain Boundary we need a (many) calculation(s) object(s)
# these should have some tag i.e. DFT-VASP-PBE. So clearly a one to many
# relationship.
# One GrainBoundary table with grain_boundary id (unique tag), name, raw atom structure
# coincident site lattice etc, and a bunch of pointers to a CalculationTable
# Each CalculationTable will have its own unique tag (DFT-VASP-PBE) along with the properties
# calculate: total energy, forces, atoms, magnetic moments.
#

grain_boundaries = {}
calculations      = {}


grain_boundaries['0000000000'] = {'title': 'Ideal Crystal ',  'gb_id':'0000000000'}
grain_boundaries['1107053111'] = {'title': 'Sigma(3)  (111)', 'gb_id':'1107053111'}
grain_boundaries['1105048332'] = {'title': 'Sigma(11) (332)', 'gb_id':'1105048332'}
grain_boundaries['1106000112'] = {'title': 'Sigma(3) (112)',  'gb_id':'1106000112'}

#
# Table energies should be populated in eV:
# Each calculation should have the atoms object attached to it.
#
calculations['0000000000'] = {'VASP-DFT-PBE': {'E0':-8.23807, 'DFT-mag': 2.2238, 'nat':1}, 'IP-EAM-MISH':{'E0': -4.2701, 'nat':1}}
calculations['1107053111'] = {'VASP-DFT-PBE' : {'E0':-406.154623782, 'nat':96, 'A': 27.7436434255}}
calculations['1105048332'] = {'IP-EAM-MISH' : {'E0':-382.847802363, 'nat':90, 'A':18.7825353894 }}
calculations['1106000112'] = {'IP-EAM-MISH' : {'E0':-196.171, 'nat':46, 'A': 9.80885920049}}
# def pretty_index(index):
# string_regex = '([0-9])(4[0-9])(3[0-9])(3[0-9])'
# a,b,c = title.findall(string_regex)
# index = 'Sigma ({0}) {1} degrees Orientation axis: [{2}] Miller plane: ({3}) '.format(a,b,c)
@app.route('/')
def home_page():
	return render_template('imeall.html', grain_boundaries=grain_boundaries.values())

@app.route('/grains/<gb_id>')
def grain_boundary(gb_id):
	gb = grain_boundaries[gb_id]
	calc = calculations[gb_id]
	print calc
	return render_template('grain_boundary.html', grain=gb, calculations=calc)
