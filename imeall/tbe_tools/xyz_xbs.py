import os
import sys
from ase import io

#color can be
# - A color is specified either as a number between 0 and 1 (gray value),
#   three numbers between 0 and 1 (red, green, blue values or RGB),
#   or as a color name from the file /usr/lib/X11/rgb.txt (or similar).

xbs_file = open("new_xbs.bs",'w')
xbs_str="atom      {}       {:.3f}      {:.3f}      {:.3f}"
#            spec     Name    Radius  Colour
spec_strs= ["spec      Fe     0.450   0.4",
            "spec      C      0.450   0.7",
            "spec      H      0.200   0.0"]
#            bonds     name 1  name 2 min-length  max-length  radius  color
bond_strs =["bonds     Fe       Fe      0.000       2.6     0.06       1.0",
            "bonds     C        Fe      0.000       2.6     0.09       0.8",
            "bonds     C        H       0.000       2.1     0.04       0.8",
            "bonds     Fe       H       0.000       2.0     0.04       1.0"]
#various parameters that can be controlled on the command line.
param_str = "inc 1"

#read xyzfile from sys.argv
ats = io.read(sys.argv[1],index="1")

print >> xbs_file, "*FeH system Migrating Fe is Labeled C"
for symbol, pos in zip(ats.get_chemical_symbols(), ats.get_positions()):
    print >> xbs_file, xbs_str.format(symbol,pos[0],pos[1],pos[2])
print >> xbs_file,""
for spec_str in spec_strs:
    print >> xbs_file, spec_str
print >> xbs_file,""
for bond_str in bond_strs:
    print >> xbs_file, bond_str
print >> xbs_file,""
print >> xbs_file, param_str
xbs_file.close()

