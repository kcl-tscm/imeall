import os
try:
  from flask import Flask, render_template
  NO_FLASK = False
except:
  NO_FLASK = True
  print 'No Flask Server Available'

def file_extension(filepath):
  endswith = filepath.split('.')[-1]
  return endswith

if not NO_FLASK:
  app             = Flask(__name__, instance_relative_config=True)
  app.config.from_pyfile('config.py')
  app.jinja_env.filters['file_extension'] = file_extension

  from imeall import views
