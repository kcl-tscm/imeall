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

The home page contains a list of links to the data view for each grain boundary
in the database. Each grain boundary page has various attributes that can be read
off and visualized.
An embedded snapshot of the grain boundary structure can be visualized,
and a 2d schematic of the grain boundary is also generated.
The database stores atomic scale properties according to calculation type
( e.g. tight binding models, EAM, DFT).

Tables of key properties: total energy, grain boundary energy,
elastic constants, Young's modulus, poisson ratio are displayed for
each calculation type.

In terms of the database structure sub-folders of each grain boundary object can
be extended to include sub calculations. Of particular relevance is a sub folder of
Nx3 atomic displacements which allows for the creation of an associated force constant
matrix at the DFT level, and defect subfolders. The design of this scheme is to make
the database suitably structured for training GAP models.

