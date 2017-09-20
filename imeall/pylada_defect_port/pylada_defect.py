from ase.calculators.vasp import Vasp
from ase.io.vasp import write_vasp

from quippy import Atoms, set_fortran_indexing
from ase import Atoms as aseAtoms
from ase.neighborlist import NeighborList
import spglib
from glob import iglob
from itertools import *
from scipy.spatial import distance
from collections import Counter
from copy import deepcopy
from quantities import eV
import numpy as np
import sys
import math
import tess
import os

set_fortran_indexing(False)

def write_interstitials(ats, ttol=0.5):
  vasp_args=dict(xc='PBE', amix=0.01, amin=0.001, bmix=0.001, amix_mag=0.01, bmix_mag=0.001,
                 kpts=[3, 3, 3], kpar=6, lreal='auto', ibrion=1, nsw=50, nelmdl=-15, ispin=2, 
                 prec='High', nelm=100, algo='VeryFast', npar=24, lplane=False, 
                 lwave=False, lcharg=False, istart=0,
                 voskown=0, ismear=1, sigma=0.1, isym=2) 
  vasp = Vasp(**vasp_args)
  vasp.initialize(ats)
  write_vasp('POSCAR', vasp.atoms_sorted, symbol_count=vasp.symbol_count, vasp5=True)

class Voronoi(object):
  """ :class:`Voronoi` is a wrapper around tess and pylada_defects methods 
  for computing and working with Voronoi tesselations.
  """
  def __init__(self, limits = ((-50,-50,-50), (50,50,50)), periodic=(False,False,False)):
    self.limits = limits
    self.periodic = periodic

  def compute_voronoi(self, points):
      """ Function to return the tess container having information of Voronoi cells 
      for given points in 3D space. 

      Args:
        points: numpy array of points coordinates.                        

      Returns:
        :tess:class:`Container`: tess container of Voronoi cells.
      """
  
      P = np.array(points)
      cntr = tess.Container(P, self.limits, self.periodic)
  
      return cntr
  
  ##########################################################
  def calculate_midpoint(self, p1, p2):
      """Calculate the midpoint given two points.
      Args:    
          p1, p2 = numpy array of point coordinates

      Returns: 
        midpoint numpy array of midpoint coordinates
      """
  
      return ((p1[0]+p2[0])/2.0, (p1[1]+p2[1])/2.0, (p1[2]+p2[2])/2.0)
  
  ##########################################################
  def calculate_polygon_centroid(self, poly_pts):
      """ Function to calculate the centroid of non-self-intersecting polygon.

      Args:
          pts: numpy array of coordinates of vertices of polygon

      Returns:
          centroid numpy array of centroid coordinates
      """
  
      P = np.array(poly_pts)
      C = np.mean(P, axis=0)
  
      return C
  
  ##########################################################
  def neighbor_list(self, list_):
      """Generator to form unique neighboring pairs along the polygon perimeter.

      Args:
        list: list of indices of Voronoi vertices forming the perimeter.

      Returns:
        list: tuple of neighboring pairs as tuples
      """
  
      i = 0
      while i + 1 < len(list_):
          yield (list_[i], list_[i+1])
          i += 1
      else:
          yield (list_[i], list_[0])
  
  ##########################################################
  def get_vertices(self, site_num, cntr):
      """Function that returns vertices of the Voronoi associated with given site.

      Args:
        site_num(int): number for the lattice site of interest
        cntr(:tess:class:`Container`): Container having information of Voronoi cells
  
      Returns:
        numpy array of Voronoi vertices coordinates
      """
  
      list_voronoi_vertices = cntr[site_num].vertices()
      V = list_voronoi_vertices
  
      # convert the list to numpy array
      V = np.asarray(V)
  
      return V
  
  ##########################################################
  def get_edgecenter(self, site_num, cntr):
      """ Function that returns vertices unique edge centers of 
      the Voronoi associated with specific lattice site.
  
      Args:
        site_num(int) for the lattice site of interest
        cntr(:tess:class:`Container`) having information of Voronoi cells
  
      Returns:
        numpy array of Voronoi edge center coordinates
      """
  
      list_face_vertices_indices = cntr[site_num].face_vertices()
  
      V_vertices = self.get_vertices(site_num, cntr)
  
      all_midpoint = []
  
      for face in list_face_vertices_indices:
          for(x,y) in self.neighbor_list(face):
              midpoint = self.calculate_midpoint(V_vertices[x], V_vertices[y])
              all_midpoint.append(midpoint)
  
      #using set so to choose only unique edge centers
      S = set(all_midpoint)
  
      #converting set to list
      Ec = list(S)
  
      #converting list to numpy array
      Ec = np.asarray(Ec)
  
      return Ec
  
  ##########################################################
  def get_facecentroid(self, site_num, cntr):
      """Function the returns vertices of face centers of the Voronoi associated with 
      specific lattice site.
  
      Args:
        site_num(int): number for the lattice site of interest.
        cntr(:tess:class:`Container`) having information of Voronoi cells.
  
      Returns:
        numpy array of Voronoi face center coordinates
      """
  
      list_face_vertices_indices = cntr[site_num].face_vertices()
  
      V_vertices = self.get_vertices(site_num, cntr)
  
      list_face_centroid = []
  
      for face in list_face_vertices_indices:
          l = []
          for j in face:
              vv  = V_vertices[j]
              l.append(vv)
          l = np.asarray(l)
          pc = self.calculate_polygon_centroid(l)
          list_face_centroid.append(pc.tolist())
  
      Fc = list_face_centroid
  
      # converting list to numpy array                                                               
      Fc = np.asarray(Fc)
  
      return Fc


def get_unique_wyckoff(ats_in):
    """ Function to find unique wyckoff sites in the primitive cell
    
    Args:
      ats_in(:ase:class:`Atoms`): Atoms object of interface.

    Returns
      unique list of lattice sites in primitive cell [[position,"chem_abbrev"],...]
    """
    # compute space group for the given primitive cell using spglib
    ats = ats_in.copy()
    #sym = spglib.get_symmetry(ats, 0.1)
    sym = spglib.get_symmetry(ats)

    # compute inverce cell of the primitive cell
    inverse_cell = np.linalg.inv(ats.cell)
    
    dummy_list = []
    wyckoff_list = []
    #args = np.argsort(ats.positions[:, 0])
    #ats = ats[args] 
    for i, at in enumerate(ats):
        print i,'/',len(ats)
        a = at.position
        a3 = get_pos_in_prim_cell(ats.copy(), a)
        frac_a = np.dot(inverse_cell, a3)
        symm_list = []
        for j in range(len(sym['rotations'])):
            # rotation matrix from sym
            R = sym['rotations'][j]
            # translation vector from sym
            Tt = sym['translations'][j]
            frac_symm_a = np.dot(R, frac_a)+Tt
            symm_a = np.dot(ats.cell, frac_symm_a)
            symm_list.append(symm_a)
            symm_a2 = get_pos_in_prim_cell(ats.copy(), symm_a)
            symm_list.append(symm_a2)
        # loop to find symmetrical equivalent positions
        for k in range(i+1, len(ats)):
            b = ats[k].position
            b2 = get_pos_in_prim_cell(ats, b)
            # check distance between positions in pos and symm_list
            if any(distance.euclidean(b, c) < 0.1 for c in symm_list):
                ats[k].position = ats[i].position
        dummy_list = a.tolist()
        dummy_list.append(ats[i].symbol)
        wyckoff_list.append(dummy_list)
    
    ### getting unique array of positions
    unique_list = [list(t) for t in set(map(tuple, wyckoff_list))]
    
    return unique_list

def get_pos_in_prim_cell(ats, a):
    """ Function to to map positions onto the primitive cell

    Parameters
        prim = pylada primitive cell
        a = cartesian coordinates of position

    Returns
        a2 = cartesian coordinates, such that fractional coordination = [0,1)
    """
    
    a1 = np.array(a)
    inv_cell = np.linalg.inv(ats.cell)

    frac_a1 = np.dot(inv_cell, a1)
    for m in range(len(frac_a1)):
        if frac_a1[m] < 0.: frac_a1[m] = frac_a1[m] + 1.
        if frac_a1[m] >= 0.999: frac_a1[m] = frac_a1[m] - 1.
    a2 = np.dot(ats.cell, frac_a1)

    return a2

def get_ints_in_prim_cell(ats, positions, interstitial_type=['B','C','N']):
    """ Function to to map positions into the primitive cell.
    
    Args:
      prim(:ase:class:`Atoms`) = pylada primitive cell
      positions(list) = list of sites (as list ['element',[x,y,z]])
      interstitial_type(list) = if type is in list, ['B' 'C' 'N'] it will
        be append to the list of sites in the prim cell

    Returns
        int_pos_list1 = unique list of sites within primitive cell
    """
   #generate copys of Atoms object and positions list. 
    prim1 = ats.copy()
    ints1 = deepcopy(positions)
    inverse_cell1 = np.linalg.inv(prim1.cell)
    int_pos_list1 = []
    
    for i in range(len(ints1)):
        a1 = np.array(ints1[i][1])
        frac_a1 = np.dot(inverse_cell1, a1)
        for m1 in range(len(frac_a1)):
            if frac_a1[m1] < 0.: frac_a1[m1] = frac_a1[m1] + 1.
            if frac_a1[m1] >= 0.999: frac_a1[m1] = frac_a1[m1] - 1.
        a2 = np.dot(prim1.cell, frac_a1)
        int_pos_list1.append([ints1[i][0], a2])
    
    return int_pos_list1

def get_interstitials(ats, ttol=0.5):
    """Function to return unique interstitial sites in the given structure.
    Args:
      Atoms(:ase:class:`Atoms`): atoms object.
      ttol(float): tolerance on distance between interstitials to be considered unique.

    Returns:
      list of unique intestital sites (cartesian (x,y,z) as list) in the given structure.
    """

#    s = deepcopy(structure)
#    prim = primitive(s)
    spg = spglib.get_spacegroup(ats, 0.1)

    ### Step 1: get unique sites in the primitive of the given structure
    uniq_sites = get_unique_wyckoff(ats)

    ### Step 2: get all interstitial sites from Voronoi method
    ints2 = get_all_interstitials(ats, uniq_sites)

    ### get interstital sites within primitive cell
    ints_prim = get_ints_in_prim_cell(ats, ints2)

    ### Step 3: get unique interstitials after symmetry analysis
    ints3 = get_unique_ints(ats, ints_prim, ttol=ttol)
        
    return ints3

def get_all_interstitials(ats, unique_atoms):
    """Function to return list of all interstitial sites using Voronoi.py

    Args:
      ats(:ase:class:`Atoms`):Atoms object of structure.
      positions(numpy array): Lattice sites to compute interstitials for.

    Returns:
      list of list, with inner list containing ['Interstitial Type', [x,y,z]]
      Interstitial type can be of type 'B', 'C', 'N' corresponding to;
        'B' Voronoi vertices.
        'C' face centers.
        'N' Edge centers.
    """

    ints_list = []
    #generate neighbours list.
    for site_num in unique_atoms:
        nl = NeighborList([2.0]*len(ats), bothways=False, self_interaction=True)
        nl.update(ats)
        site = np.array([ats[site_num].position[0], ats[site_num].position[1], ats[site_num].position[2]])
        points = []
        indices, offsets = nl.get_neighbors(site_num)
        print 'site_num:', site_num
        for i, offset in zip(indices, offsets):
            if offset[2] == 0.0:
              print 'neighb:', i,ats.positions[i], offset, np.dot(offset, ats.get_cell())
              points.append(ats.positions[i] + np.dot(offset, ats.get_cell()))
        ### converting list to numpy array
        print "neighbors", len(points)
        points = np.asarray(points)

        voronoi = Voronoi()

        ### using tess object cntr to compute voronoi
        cntr = voronoi.compute_voronoi(points)

        ### Voronoi vertices
        ### the first position in points is the site, therefore '0'
        v = voronoi.get_vertices(0, cntr)

        for i in range(len(v)):
            ints_list.append(['B', v[i].tolist()])

        ### Voronoi face centers
        f = voronoi.get_facecentroid(0, cntr)

        for j in range(len(f)):
            ints_list.append(['C', f[j].tolist()])

        ### Voronoi edge centers
        e = voronoi.get_edgecenter(0, cntr)

        for k in range(len(e)):
            ints_list.append(['N', e[k].tolist()])

    ### return list of list ['Atom type', [x,y,z]]
    return ints_list


def get_unique_ints(ats, int_pos, ttol=0.5):
    """ Function to find unique interstitial sites in the primitive cell.

    Args:
      ats(:ase:class:`Class`)  Atoms object
      int_pos(list): of interstitial positions in primitive cell
      ttol(float): tolerance, default = 0.5

    Returns:
      list of unique intestital sites (cartesian (x,y,z) as list) in primitive cell
    """

    pos2 = deepcopy(int_pos)
    sym2 = spglib.get_symmetry(ats)

    inverse_cell2 = np.linalg.inv(ats.cell)
    int_list2 = []

    for i in range(len(pos2)):
    #   numpy 1D array into column vector with shape listed as (3,)
        print '{}/{}'.format(i, len(pos2))
        a4 = pos2[i]
        frac_a4 = np.dot(inverse_cell2, a4)
        symm_list2 = []
        for j in range(len(sym2['rotations'])):
            R = sym2['rotations'][j]
            Tt = sym2['translations'][j]
            frac_symm_a4 = np.dot(R, frac_a4) + Tt
            symm_a4 = np.dot(ats.cell, frac_symm_a4)
            symm_list2.append(symm_a4)
            symm_a5 = get_pos_in_prim_cell(ats.copy(), symm_a4)
            symm_list2.append(symm_a5)
        # loop to find symmetrical equivalent positions    
        for k in range(i+1, len(pos2)):
            b2 = pos2[k]
            # check distance between positions in pos and symm_list
            if any(distance.euclidean(b2,c2) < ttol for c2 in symm_list2):
                pos2[k] = pos2[i]
        int_list2.append(pos2[i])

    ### getting the unique list of interstitials
    unique_int_list = [list(t) for t in set(map(tuple, int_list2))]

    return unique_int_list

