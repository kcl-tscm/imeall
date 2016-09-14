import os
try:
  from flask import Flask, render_template
except:
  print 'No Flask Server Available'

if __name__=='__main__':
  app             = Flask(__name__, instance_relative_config=True)
  app.config.from_pyfile('config.py')
  from imeall import views
