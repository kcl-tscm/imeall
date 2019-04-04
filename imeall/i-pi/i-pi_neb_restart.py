import os
import sys
import glob
import numpy as np
from ase import io
from ase.io.extxyz import write_xyz

def write_restart(images):
    conv = 1.0000
    alat = 14.3499
    header_string = "# CELL(abcABC):   {alat}    {alat}    {alat}    90.00000    90.00000    90.00000  Step:           0  Bead: {bead} positions{{angstrom}}  cell{{angstrom}}"
    with open('neb_restart.xyz','w') as f:
        for bead, ats in enumerate(images):
            if bead == 0:
                ats.write('posi.xyz')
            if bead == len(images)-1:
                ats.write('posf.xyz')

            print >> f, len(ats)
            print >> f, header_string.format(alat=alat, bead=bead+1)
            for symbol, pos in zip(ats.get_chemical_symbols(), ats.get_positions()):
                #to fix small discrepancy in lattice constants.
                conv_pos = conv*pos
                print >> f, "       {}  {} {} {}".format(symbol, conv_pos[0], conv_pos[1], conv_pos[2])

prefix = sys.argv[1]
positions = glob.glob("{prefix}.pos_*xyz".format(prefix=prefix))
func = lambda x: int(x.split("_")[1].split(".")[0])
positions.sort(key = func)

if os.path.isfile("neb_restart.xyz"):
    os.remove("neb_restart.xyz")

images = []
for pos_file in positions:
    ats = io.read(pos_file,index="-1")
    images.append(ats)
write_restart(images)
