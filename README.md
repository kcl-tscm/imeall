`Imeall` is a database framework for the calculation of 
the atomistic properties of grain boundaries.
The grain boundaries are generated and labeled according to 
a uniform process  indicated by the following notation:

```
  Theta[XXX](YYY)
```

sigma = (number of coincident sites in an elementary cell)/
(total number of all lattice sites in an elementary cell)

Theta is the angle of misorientation about the vector [XXX]
and (YYY) determines the plane of the interface between
the misoriented Grain A and Grain B. In this way it is possible
to uniquely characterize all grain boundaries. 

An additional set of three degrees of freedom is then used to extend 
the above definition. These degrees of freedom correspond to a 
rigid body translation in the two dimensional plain of the 
grain boundary one relative to the other, and a volume relaxation
orthogonal to the grainboundary plane. These atomic scale degrees
of freedom compose subgrains of the canonical parent grain.
Finally, we assign the ideallized crystalline bulk structure the notation:

```
  0000000000
```

Computationally each grain boundary so identified can be considered an
object with all its physical properties as attributes (to be described below).
The structures are stored locally in the database as [Extended
XYZ](https://libatoms.github.io/QUIP/io.html#module-ase.io.extxyz) files.

## Dependencies
  - [Flask](http://flask.pocoo.org/)
  - [ASE](https://wiki.fysik.dtu.dk/ase/)
  - [QUIP](https://libatoms.github.io/QUIP/quippy.html)
  - [OVITO](https://www.ovito.org)

## Using Imeall
0. To get the development branch of Imeall:
```
  git clone https://github.com/Montmorency/imeall.git
```

Set the environment variable (preferably in your .bashrc or .bash_profile)
for your ovito executable and the directory containing any .xml files required
for QUIP potential calculators.
```
  export POTDIR=~/imeall/imeall/potentials
  export OVITO=/Applications/ovito
```

1. Initialize the imeall app:
```
  python runserver.py
```
This will launch the web app hosted locally and can be accessed through
your browser at http://127.0.0.1:5000/.

The home page contains a list of links for each material
in the database. Each link leads to a view of different orientation axes vectors for 
which grain boundary properties have been computed. 

For instance /alphaFe/110 contains all computed 110 grain boundaries for iron. The 000 
orientation axes implies a perfect crystal, different cleavage planes 
and singular defects and vacancies in a perfect crystal. 

Within each orientation axes is a list of full grain boundaries. 
Each of these contains a snap shot of the unrelaxed grain boundary structure
(the canonical boundary), the coincident site lattice, and a
list of sub-directories containing the type of calculations performed on the
grain boundary. The subdirectories (calc_types) store atomic scale 
properties according to calculation type. Examples of calculation types include:
tight binding models, EAMs, DFT. If multiple "flavours" of each calculation
type are to be used the directories can take more specific 
names: EAM_Mishin, DFT_PW, DFT_VASP etc.
Within each calc_type of directories is the pattern 
gbid_'suffix'. The suffix denotes the modification to the canonical grain boundary
structure. This could constitute a pattern like ```_d2.0``` which denotes:
delete one of a pair of atoms from the original grain 
boundary structure with nearest neighbour distance less than 2 A.
Another suffix pattern might be displace the canonical grain along the in-plane
lattice vectors by a distance 0.25*lattice vector:
_tv0.25bxv0.25. Suffices can be combined ```_tv0.25bxv0.25_d2.0```.
A key for some common suffix patterns is given here:
```
  _v2bxv6        : supercell 2 times along v 6 times along bxv.
  _tv0.13bxv0.25 : displace one grain along the 0.13*v along v and 0.13*bxv along bxv.
  _d2.0z         : delete one of pair of atoms within distance 2.0 A of each other.
  _h2.0n1        : 1 additional hydrogen atom in system 2.0 A from a grainboundary.
```

Systematic use of these naming conventions ensures that at the analysis stage
grain boundary properties can be retrieved, compared, and combined easily.

The canonical grain contains a file: 
```
  gb.json
```
subsequent directories each contain a file:
```
  subgb.json
```
This json file contains information about relevant grain boundary properties
in the usual key:value format. Some common keys are: 
```
  E_gb_init : Unrelaxed Grain Boundary Energy.
  E_gb      : Relaxed Grain Boundary Energy.
  n_at      : Number of atoms in the grain boundary.
  A         : Grain Boundary area in the unit cell.
```
Upon subsequent analysis new key:values pairs can be added to the 
json dictionaries as required.

