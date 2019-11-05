import os
import glob
import numpy as np

from ase.io import read
from ase.io.xyz import write_xyz
from argparse import ArgumentParser
from quippy import Atoms, set_fortran_indexing

set_fortran_indexing(False)

#to select qm atoms from the mm_cell.
def pull_file(args, dir_name=None):
    if dir_name == None:
        ats = read(args.input, index=":")
    else:
        ats = read(os.path.join(dir_name, args.input), index=":")

    if not args.full_qm:
        init_ats = Atoms("crack_sim.xyz")
        crack_pos = init_ats.info["CrackPos"]

        qm_radius = 3.0
        buff = 8.0
        mm_pos = init_ats.get_positions()
        x, y, z = init_ats.positions.T
#radius1 = np.sqrt((x - crack_pos[0])**2 + (y-crack_pos[1])**2 + (z-crack_pos[2])**2)
        radius1 = np.sqrt((x - crack_pos[0])**2 + (y-crack_pos[1])**2)
        qm_region_mask = (radius1 < qm_radius)
    
        qm_pos = init_ats.get_positions()[qm_region_mask]
        qm_com = qm_pos.mean(axis=0) # just to avoid pbc errors..

        cell_center = np.diag(ats[0].get_cell())/2.0
        print cell_center
        print 'Full Cell CrackPos: ', crack_pos
        crack_pos = crack_pos - qm_com + cell_center
        print 'Shift CrackPos', crack_pos
        print 'QM Radius', qm_radius
        print len(ats)
        for at in ats[args.start:args.end]:
       #mark quantum atoms
            x,y,z = at.positions.T
            radius1 = np.sqrt((x - crack_pos[0])**2 + (y-crack_pos[1])**2 + (z-crack_pos[2])**2)
            qm_region_mask = np.array(map(int, (radius1 <= qm_radius)))
            print sum(qm_region_mask)
            at.new_array("qm_atoms", qm_region_mask, dtype=float)
            at.set_array("qm_atoms", qm_region_mask)
            write_xyz(args.output, at, append=True)
    else:
        for at in ats[args.start:args.end]:
            write_xyz(args.output, at, append=True)
        

parser = ArgumentParser()
parser.add_argument("--start", "-s", type=int, default=0, help="for multiple force iterations in output file start range.")
parser.add_argument("--end", "-e", type=int, default=1, help="for multiple force iterations in output set end number.")
parser.add_argument("--input", "-i", default="OUTCAR")
parser.add_argument("--output", "-o", default="forces.xyz")
parser.add_argument("--full_qm", "-f", action="store_true")
parser.add_argument("--pattern", "-p", help="Job directory patterns.", default=None)
parser.add_argument("--qm_radius", "-q", default=3.0, type=float)
args = parser.parse_args()

if args.pattern != None:
    j_dirs = glob.glob(args.pattern)
    for j_dir in j_dirs:
        print j_dir
        pull_file(args, dir_name=j_dir)
else:
    pull_file(args, dir_name=None)
