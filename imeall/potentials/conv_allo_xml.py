import os
import sys
from ase.calculators.eam import EAM
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as spline


class EAMXML(object):
  def __init__(self, input_file='undefined', rmin = 0.0, rcut=5.3, ds=1):
    self.input_file  = input_file
    self.output_file = 'quip.xml'
    self.rmin        = rmin
    self.rcut        = rcut
    self.N           = 26
    self.ds          = ds
    self.drho        = 0.029

  def plot_xml(self, xml_file_name):
    with open(xml_file_name, 'r') as f:
      root = ET.parse(f)
    V = [(x.attrib['r'], x.attrib['y'] ) for x in root[1][0]]
    print V

  def alloy_to_xml(self, n_spec=1, F=[], V=[], rho=[]):
  #   Converts EAM potential to quippy friendly xml. n_spec number of species.
  #    F embedding terms for each atom, rho, density terms, V pair-potentials.
    at_types = [26,1]
    at_pairs = [[26,26],[1,26],[1,1]]
#Write initial Densities and Embedding Terms for each species.
    i_spec = 0
    dat_F   = len([i for i,f in enumerate(F[:,i_spec]) if i%self.ds==0])
    dat_rho = len([i for i,f in enumerate(V[:,i_spec]) if i%self.ds==0])
    dat_V   = len([i for i,f in enumerate(rho[:,i_spec]) if i%self.ds==0])
    qeam_file = open('quip.xml','w')

#Write EMBEDDING and DENSITY FUNCTIONS:
    print >> qeam_file, "<EAM_ErcolAd_params n_types=\"{3}\" n_spline_V=\"{2}\" n_spline_rho=\"{0}\" n_spline_F=\"{1}\">".format(dat_rho, dat_F, dat_V, n_spec)
    for i_spec in range(n_spec):
      print i_spec
      nsteps  = 10000
      r       = [float(i)*self.rcut/float(nsteps) for i in range(nsteps)]
      fr      = [self.drho*float(i) for i in range(0, nsteps)]
      print >> qeam_file, "<per_type_data atomic_num=\"{0}\" type=\"{1}\">".format(at_types[i_spec], i_spec+1)
#DENSITY TERM
      print >> qeam_file, "<spline_rho>"
      for i, d in enumerate(zip(r, rho[:,i_spec])): 
        if d[0] >= self.rmin and i%self.ds==0:
          print >> qeam_file, '\t\t\t<point r=\"{0:.7f}\" y=\" {2:.7f} \" b=\" {1} \" c=\" {1} \" d=\" {1} \"/>'.format(d[0], '0.0000000', d[1])
      print >> qeam_file, "</spline_rho>"
#EMBEDDING TERM
      print >> qeam_file, "<spline_F>"
      for i, d in enumerate(zip(fr, F[:,i_spec])): 
        if i%self.ds==0:
          print >> qeam_file, '\t\t\t<point r=\"{0}\" y=\" {2} \" b=\" {1} \" c=\" {1} \" d=\" {1} \"/>'.format(d[0], '0.0000000', d[1])
      print >> qeam_file, '</spline_F>'
      print >> qeam_file, "</per_type_data>"

#Now write PAIR POTENTIALS: 
    for i_pair in range(n_pairs):
      type_i = at_pairs[i_pair][0]
      type_j = at_pairs[i_pair][1]
      print >> qeam_file, "<per_pair_data atomic_num_i=\"{2}\" atomic_num_j=\"{3}\" r_min=\"{0}\" r_cut=\"{1}\">".format(self.rmin, self.rcut, type_i, type_j)
      print >> qeam_file, "<spline_V>"
      phi = spline(r, V[:,i_pair], k=3)
      for i, d in enumerate(r): 
        if d>=self.rmin and i%self.ds==0:
          print >> qeam_file, '\t\t\t<point r=\"{0}\" y=\" {2} \" b=\" {1} \" c=\" {1} \" d=\" {1} \"/>'.format(d, '0.0000000', 0.5*V[i,i_pair])
      print >> qeam_file, "</spline_V>"
      print >> qeam_file, "</per_pair_data>"

    at_pairs = [[26,1],[1,26]]
#Now write PAIR Density Splines: 
    #for i_pair in range(2):
     # type_i = at_pairs[i_pair][0]
     # type_j = at_pairs[i_pair][1]
     # print >> qeam_file, "<per_pair_data atomic_num_i=\"{2}\" atomic_num_j=\"{3}\" r_min=\"{0}\" r_cut=\"{1}\">".format(self.rmin, self.rcut, type_i, type_j)
     # print >> qeam_file, "<spline_rho>"
     # for i, d in enumerate(zip(r, rho[:,2+i_pair])): 
     #   if d[0] >= self.rmin and i%self.ds==0:
     #     print >> qeam_file, '\t\t\t<point r=\"{0:.7f}\" y=\" {2:.7f} \" b=\" {1} \" c=\" {1} \" d=\" {1} \"/>'.format(d[0], '0.0000000', d[1])
     # print >> qeam_file, "</spline_rho>"
     # print >> qeam_file, "</per_pair_data>"
    print >> qeam_file, "</EAM_ErcolAd_params>"

if __name__ == '__main__':
#If we convert the file into the Fe.alloy format then ase can read it directly.
  eam_parse = EAMXML(rmin=0.0, rcut=5.3, ds=10)
  PLOT_POT = False
  CONV_TO_XML = True
  n_spec   = 2
  n_pairs  = n_spec*(n_spec+1)/2
  print n_spec, n_pairs
  n_points = 10000
#Atom Embedding Term
  F   = np.zeros([n_points, n_spec])
#Atom Pair Potential
  V   = np.zeros([n_points, n_pairs])
#Atom Pair Potential
  rho = np.zeros([n_points, 4])
#Uncomment to read in Dudarev potential:
  if CONV_TO_XML:
    F[:,0]   = map(float, open('./FeHPotB/Fe_F.dat','r').read().split())
    V[:,0]   = map(float, open('./FeHPotB/Fe_Fe_V.dat','r').read().split())
    rho[:,0] = map(float, open('./FeHPotB/Fe_Fe_rho.dat','r').read().split())
    F[:,1]   = map(float, open('./FeHPotB/H_F.dat','r').read().split())
    V[:,1]   = map(float, open('./FeHPotB/H_H_V.dat','r').read().split())
    rho[:,1] = map(float, open('./FeHPotB/H_H_rho.dat','r').read().split())
    V[:,2]   = map(float, open('./FeHPotB/Fe_H_V.dat','r').read().split())
    rho[:,2] = map(float, open('./FeHPotB/Fe_H_rho.dat','r').read().split())
    rho[:,3] = map(float, open('./FeHPotB/H_Fe_rho.dat','r').read().split())
    eam_parse.alloy_to_xml(n_spec=2, F=F, V=V, rho=rho)

  if PLOT_POT:
    with open('quip.xml', 'r') as f:
      tree = ET.parse(f)
    root = tree.getroot()
    R = [float(x.attrib['y']) for x in root[0][0]]
    F = [float(x.attrib['y']) for x in root[0][1]]
    V = [float(x.attrib['y']) for x in root[4][0]]
    xpts = []
    for x in root.findall('./per_type_data/spline_rho/'):
      xpts.append(x.attrib['r'])
    xpts = xpts[:len(R)]
    rho = [float(x.attrib['r']) for x in root[0][1]]
    print 'length of vectors', len(xpts), len(rho)
    fig, ax = plt.subplots(3,1)
    ax[0].set_ylim([-1,2])
    ax[0].set_ylabel('Density')
    ax[0].plot(xpts, R)
    ax[1].set_ylim([-5,50])
    ax[1].set_ylabel('V')
    ax[1].plot(xpts,V)
    ax[2].set_ylim([-50,20])
    ax[2].set_ylabel('Embedding Term')
    ax[2].plot(rho, F)
  plt.show()
