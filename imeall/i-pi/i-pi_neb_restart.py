import os
import sys
import glob
import numpy as np
from ase import io
from ase.io.extxyz import write_xyz


prefix = sys.argv[1]
positions = glob.glob("{prefix}.pos_*xyz".format(prefix=prefix))
positions.sort(key = func)
func = lambda x: int(x.split("_")[1].split(".")[0])

if os.path.isfile("neb_restart.xyz"):
    os.remove("neb_restart.xyz")

for pos_file in positions:
    ats = io.read(pos_file,index="-1")
    write_xyz("neb_restart.xyz",ats, append=True)

