import os
import sys
from ase.calculators.eam import EAM
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as spline


class EAMXML(object):
  def __init__(self, input_file='undefined', rmin = 1.0, rcut=5.3, ds=1):
    self.input_file  = input_file
    self.output_file = 'quip.xml'
    self.rmin        = rmin
    self.rcut        = rcut
    self.N           = 26
    self.ds          = ds
    self.drho        = 0.0023

  def plot_xml(self, xml_file_name):
    with open(xml_file_name, 'r') as f:
      root = ET.parse(f)
    V = [(x.attrib['r'], x.attrib['y'] ) for x in root[1][0]]
    print V

  def alloy_to_xml(self, FVR=None):
    if FVR == None:
      print 'INPUT_FILE', self.input_file
      alloy     = EAM(potential=self.input_file)
      qeam_file = open('quip.xml','w')
      dat_rho = len([x for i, x in enumerate(zip(alloy.r, alloy.density_data[0])) if x[0] >= self.rmin and i%self.ds==0])
      dat_F   = len([x for i, x in enumerate(alloy.rho) if i%self.ds==0])
      dat_V   = len([d for i, d in enumerate(alloy.r) if d >= self.rmin and i%self.ds==0])
      r   = alloy.r
      fr  = alloy.rho
      rho = alloy.density_data[0]
      F   = alloy.embedded_data[0]
      phi = alloy.phi[0,0]
    else:
      nsteps  = 10000
      dat_F   = len([i for i,f in enumerate(FVR[0]) if i%self.ds==0])
      dat_rho = len([i for i,f in enumerate(FVR[1]) if i%self.ds==0])
      dat_V   = len([i for i,f in enumerate(FVR[2]) if i%self.ds==0])
      qeam_file = open('quip.xml','w')
      r       = [float(i)*self.rcut/float(nsteps) for i in range(0, nsteps)]
      fr      = [self.drho*float(i) for i in range(0,nsteps)]
      #F       = [f for i,f in enumerate(FVR[0]) if i%self.ds==0]
      #V       = [f for i,f in enumerate(FVR[1]) if i%self.ds==0]
      #rho     = [f for i,f in enumerate(FVR[2]) if i%self.ds==0]
      F = FVR[0]
      V = FVR[1]
      rho = FVR[2]
      phi = spline(r, V, k=3)

    print >> qeam_file, "<EAM_ErcolAd_params n_types=\"1\" n_spline_V=\"{2}\" n_spline_rho=\"{0}\" n_spline_F=\"{1}\">".format(dat_rho, dat_F, dat_V)
    print >> qeam_file, "<per_type_data atomic_num=\"26\" type=\"1\">"
    print >> qeam_file, "<spline_rho>"
    for i, d in enumerate(zip(r,rho)): 
      if d[0] >= self.rmin and i%self.ds==0:
        print >> qeam_file, '\t\t\t<point r=\"{0:.7f}\" y=\" {2:.7f} \" b=\" {1} \" c=\" {1} \" d=\" {1} \"/>'.format(d[0], '0.0000000', d[1])
    print >> qeam_file, "</spline_rho>"

    print >> qeam_file, "<spline_F>"
    for i, d in enumerate(zip(fr,F)): 
      if i%self.ds==0:
        print >> qeam_file, '\t\t\t<point r=\"{0}\" y=\" {2} \" b=\" {1} \" c=\" {1} \" d=\" {1} \"/>'.format(d[0], '0.0000000', d[1])
    print >> qeam_file, '</spline_F>'
    print >> qeam_file, "</per_type_data>"
    print >> qeam_file, "<per_pair_data atomic_num_i=\"{0}\" atomic_num_j=\"{0}\" r_min=\"{1}\" r_cut=\"{2}\">".format(self.N, self.rmin, self.rcut)
    print >> qeam_file, "<spline_V>"
    for i, d in enumerate(r): 
      if d>=self.rmin and i%self.ds==0:
        print >> qeam_file, '\t\t\t<point r=\"{0}\" y=\" {2} \" b=\" {1} \" c=\" {1} \" d=\" {1} \"/>'.format(d, '0.0000000', phi(d))
    print >> qeam_file, "</spline_import xml.etree.ElementTree as ETV>"
    print >> qeam_file, "</per_pair_data>"
    print >> qeam_file, "</EAM_ErcolAd_params>"

if __name__ == '__main__':
#If we convert the file into the Fe.alloy format then ase can read it directly.
  eam_parse = EAMXML(rmin=0.0, rcut=5.3, ds=4)
  FVR = []
#Uncomment to read in Dudarev potential:
  #F = map(float, open('./Dudarev/fort.1000','r').read().split())
  #V = map(float, open('./Dudarev/fort.1002','r').read().split())
  #R = map(float, open('./Dudarev/fort.1001','r').read().split())
  #print len(F)
  #print len(V)
  #print len(R)
  #FVR.append(F)
  #FVR.append(V)
  #FVR.append(R)
  #eam_parse.alloy_to_xml(FVR)
  with open('Fe_Dudarev.xml', 'r') as f:
    tree = ET.parse(f)
  root = tree.getroot()
  R = [float(x.attrib['y']) for x in root[0][0]]
  F = [float(x.attrib['y']) for x in root[0][1]]
  V = [float(x.attrib['y']) for x in root[1][0]]
  xpts = []
  for x in root.findall('./per_type_data/spline_rho/'):
    xpts.append(x.attrib['r'])
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
