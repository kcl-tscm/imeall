imeall is a database of atomistic properties, structure, energy, and forces,
related to grain boundaries in alpha-Fe.

The grain boundaries are generated and labeled according to 
a uniform process  indicated by the following notation:

sig(A)Theta[XXX](YYY)

A is the reciprocal density of the coincident sites
in the Coincident Site lattice i.e.:

sigma = (number of coincident sites in an elementary cell)/
(total number of all lattice sites in an elementary cell)

Theta is the angle of misorientation about the vector [XXX]
and (YYY) determines the plane of the interface between
the misoriented grain b and grain a. In this way it is possible
to uniquely characterize all the grain boundaries. 

An additional set of three degrees of freedom is then used to extend 
the above definition. These degrees of freedom correspond to a 
rigid body translation of one grain relative to the other. These
however are limited to rigid translation the result in stable 
atomic configurations of the grain boundary. Hence for the purposes
of this database the Grain boundary is determined by the first
five degrees of freedom, within each of the Grainboundaries there is
then a possibility to represent different rigid body translations.

We assign the ideallized crystalline bulk structure the notation:
\sigma(0)0[ooo](nnn)

Computationally each grain boundary so identified can be considered an
object with all its physical properties as attributes (to be described below).
The structures are stored locally in the database as ASE atoms objects.
For access and comparison the grain boundary objects are stored as
a dictionary with key 'CCCCooonnn' i.e. CCCC the angular misorientation
in degrees (4 decimal places should be sufficient) AAA, the
rotation axis ooo, and nnn the miller index of the boundary plane.

## Dependencies
  - [Flask](http://flask.pocoo.org/)
  - [ASE](https://wiki.fysik.dtu.dk/ase/)

1. Initialize the imeall app:
```
    python runserver.py
```
This will launch the web app hosted locally and can be accessed through
your browser at http://127.0.0.1:5000/.

The home page contains a list of links to the data view for each material
in the database. This leads to a view of different orientation axes vectors for 
which grain boundary properties have been computed. For instance /alphaFe/110 
contains all 110 grain boundaries for iron etc. The 000 
orientation axes implies a perfect crystal, and different cleavage planes 
and singular defects. Within each orientation axes is a list of 
full grain boundaries. Each of these should contain a snap shot
of the unrelaxed grain boundary structure, the coincident site lattice, and a
list of sub-directories containing the type of calculations performed on the
grain boundary. The subdirectories (subgrains) stores atomic scale 
properties according to calculation type (e.g. tight binding models, 
EAM, DFT, EAM_Mishin, DFT_PW, DFT_VASP) 
and within each of them a series of directories with the pattern 
gbid_'suffix'. The suffix denotes the modification to the base grain
structure. This could constitute a pattern like _d2.0. Which denotes:
delete one of a pair of atoms from the original grain 
boundary structure with nearest neighbour distance less than 2 A.
Another suffix pattern might be displace the grain along x by distance 0.25:
_tx0.25. These can be combined _d2.0_tx0.25. A key for these patterns is given
here:

_d2.0   : delete one of pair of atoms with distance < 2,0 A.
_tx2,0  : displace one grain along x 2.0 A.
_ty2.0  : displace one grain along y 2,0 A.
_vn123  : vacany created by deleting atom 123.
_h1     : 1 additional hydrogen atom in system.


## Installation
To get the development branch of Imeall:
```
  git clone https://github.com/Montmorency/imeall.git
```
