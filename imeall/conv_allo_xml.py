import os
import sys
from ase.calculators.eam import EAM
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as spline


class EAMXML(object):
  def __init__(self, input_file='./Mendelev/Mendelev.Fe.EAM.alloy', rmin = 1.0, rcut=5.3):
    self.input_file  = input_file
    self.output_file = 'eam_quip.xml'
    self.rmin        = rmin
    self.rcut        = rcut 
    self.N           = 26
    self.ds          = 4

  def plot_xml(self, xml_file_name):
    with open(xml_file_name, 'r') as f:
      root = ET.parse(f)
    V = [(x.attrib['r'], x.attrib['y'] ) for x in root[1][0]]
    print V

  def alloy_to_xml(self, FVR=None):
    if FVR == None:
      print self.input_file
      alloy     = EAM(potential=self.input_file)
      qeam_file = open('quip.xml','w')
      dat_rho = len([x for i, x in enumerate(zip(alloy.r, alloy.density_data[0])) if x[0] >= self.rmin and i%self.ds==0])
      dat_F   = len([x for i, x in enumerate(alloy.rho) if i%self.ds==0])
      dat_V   = len([d for i, d in enumerate(alloy.r) if d >= self.rmin and i%self.ds==0])
      r   = alloy.r
      fr = alloy.rho
      rho = alloy.density_data[0]
      F = alloy.embedded_data[0]
      phi = alloy.phi[0,0]
    else:
      nsteps  = 10000
      dat_F   = len(FVR[0])
      dat_rho = len(FVR[1])
      dat_V   = len(FVR[2])
      qeam_file = open('quip.xml','w')
      r       = [float(i)*self.rcut/float(nsteps) for i in range(nsteps)]
      fr      = [0.0299999*float(i) for i in range(0,10000)]
      F       = FVR[0]
      V       = FVR[1]
      rho     = FVR[2]
      phi     = spline(r, V, k=3)

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
    print >> qeam_file, "</spline_V>"
    print >> qeam_file, "</per_pair_data>"
    print >> qeam_file, "</EAM_ErcolAd_params>"

if __name__ == '__main__':
  eam_parse = EAMXML(rmin=1.0, rcut=5.30,input_file='./AcklandFeP/Fe.alloy')
  #eam_parse = EAMXML(rmin=0.0, rcut=5.30,input_file='./AcklandFeP/Fe-P.eam.alloy')
  FVR = []
  F = map(float, open('./AcklandFeH/Fe_F.dat','r').read().split())
  V = map(float, open('./AcklandFeH/Fe_Fe_V.dat','r').read().split())
  R = map(float, open('./AcklandFeH/Fe_rho.dat','r').read().split())
  print len(F)
  print len(V)
  print len(R)
  FVR.append(F)
  FVR.append(V)
  FVR.append(R)
  eam_parse.alloy_to_xml()
#  eam_parse.alloy_to_xml()
# with open('Fe_Mendelev_Untruncated.xml', 'r') as f:
# with open('iron_mish.xml', 'r') as f:
  with open('quip.xml', 'r') as f:
    tree = ET.parse(f)
  root = tree.getroot()
  R = [ float(x.attrib['y']) for x in root[0][0]]
  F = [ float(x.attrib['y']) for x in root[0][1]]
  V = [ float(x.attrib['y']) for x in root[1][0]]
  fig, ax = plt.subplots(3,1)
  ax[0].set_ylabel('Density')
  ax[0].plot(R)
  ax[1].set_xlim([2000,3000])
  ax[1].set_ylim([-50,50])
  ax[1].set_ylabel('V')
  ax[1].plot(V)
  ax[2].set_ylabel('Embedding Term')
  ax[2].plot(F)
  plt.show()
