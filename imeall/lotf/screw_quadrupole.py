import sys
import numpy as np
from ase import Atoms, Atom
from quippy import Atoms as quipAtoms
from quippy.structures import disloc_noam
from quippy import set_fortran_indexing

set_fortran_indexing(False)


alat = 2.82893

#supercell
a1 = (1./3.)*alat*np.array([-1,-1,2])
a2 = (1./2.)*alat*np.array([1,-1,0])
a3 = (1./2.)*alat*np.array([1,1,1])

#latt_1 = alat*np.array([-1,-1,2])
#latt_2 = alat*np.array([1,-1,0])
#latt_1 = alat*np.array([-1,-1,2])
#latt_2 = alat*np.array([1,-1,0])
#latt_3 = alat*np.array([1,1,1])

latt_1 = alat*np.array([1,0,0])
latt_2 = alat*np.array([0,1,0])
latt_3 = alat*np.array([0,0,1])

basis_1 = alat*np.array([0,0,0])
basis_2 = alat*np.array([0.5,0.5,0.5])
basis_3 = alat*np.array([-0.5, 0.5, -0.5])


sup_cell = 'sup_2'
core_type = "easy"

if sup_cell =='sup_1':
    n = 15.
    m = 9.
elif sup_cell=='sup_2':
    n = 21.
    m = 13.
else:
    sys.exit("not valid cell.")

if core_type == "easy":
  c1 = (n*a1) - (1.0/(3.0*m))*a3
  c2 = (n/2.)*a1 + m*a2 + (0.5)*a3 - 1.0/(6.0*m)*a3
  c3 = a3
elif core_type == "hard":
  c1 = (n*a1) + (1.0/(3.0*m))*a3
  c2 = (n/2.)*a1 + m*a2 + (0.5)*a3 + 1.0/(6.0*m)*a3
  c3 = a3

#ats = Atoms(cell=[c1,c2,c3], pbc=[True,True,True])
ats = Atoms(pbc=[True,True,True])

for x in range(-25,25):
  for y in range(-25,25):
    for z in range(-25,25):
      latt = x*latt_1 + y*latt_2 + z*latt_3
      Fe_1 = basis_1 + latt
      Fe_2 = basis_2 + latt
      Fe_3 = basis_3 + latt

      point_1 = Fe_1 
      point_2 = Fe_2 
      c3xc1 = np.cross(c3, c1)
      c3xc1 = c3xc1/(np.linalg.norm(c3xc1)**2)

      c2xc3 = np.cross(c2, c3)
      c2xc3 = c2xc3/(np.linalg.norm(c2xc3)**2)

      c1xc2 = np.cross(c1, c2)
      c1xc2 = c1xc2/(np.linalg.norm(c1xc2)**2)

      if ((np.dot(point_1, c3xc1) > 0) and (np.dot(point_1,c2xc3) > 0) and (np.dot(point_1, c1xc2) > 0) 
         and (np.dot((point_1-c3), c1xc2) <=0) and (np.dot((point_1-c2), c3xc1) <= 0) and (np.dot((point_1-c1), c2xc3) <=0)):
        ats.append(Atom('Fe', position=Fe_1))

      if ((np.dot(point_2, c3xc1) > 0) and (np.dot(point_2,c2xc3) > 0) and (np.dot(point_2, c1xc2) > 0) 
         and (np.dot((point_2-c3), c1xc2) <=0) and (np.dot((point_2-c2), c3xc1) <= 0) and (np.dot((point_2-c1), c2xc3) <=0)):
        ats.append(Atom('Fe', position=Fe_2))

if sup_cell =="sup_1":
    ats.append(Atom('Fe', position=0.5*a3) )

ats = quipAtoms(ats)
ats.set_cutoff(3.0)
ats.calc_connect()
core = 0.25*c1+ 0.5*c2 + 0.5*c3
brg_vec = -c3
disloc_l = c3
disloc_noam(ats, core, disloc_l, brg_vec)

core = 0.75*c1+ 0.5*c2 + 0.5*c3
brg_vec = +c3
disloc_l = c3
disloc_noam(ats, core, disloc_l, brg_vec)

ats.set_cell([c1,c2,c3], scale_atoms=False)
ats.write("screw.xyz")

