from ase.io import read
from ase.io.xyz import write_xyz
from argparse import ArgumentParser
from quippy import Atoms, set_fortran_indexing

import numpy as np

set_fortran_indexing(False)

parser = ArgumentParser()
parser.add_argument("--start", "-s", type=int, default=0)
parser.add_argument("--input", "-i", default="OUTCAR")
parser.add_argument("--output", "-o", default="forces.xyz")
parser.add_argument("--qm_radius", "-q", default=3.0, type=float)
args = parser.parse_args()

init_ats = Atoms("crack_sim.xyz")
crack_pos = init_ats.info["CrackPos"]

#to select qm atoms from the mm_cell.
qm_radius = 3.0
buff = 8.0
mm_pos = init_ats.get_positions()
x, y, z = init_ats.positions.T

#radius1 = np.sqrt((x - crack_pos[0])**2 + (y-crack_pos[1])**2 + (z-crack_pos[2])**2)
radius1 = np.sqrt((x - crack_pos[0])**2 + (y-crack_pos[1])**2)

qm_region_mask = (radius1 < qm_radius)
#qm_buffer_mask = (radius1 < qm_radius + buff)
qm_pos = init_ats.get_positions()[qm_region_mask]
qm_com = qm_pos.mean(axis=0) # just to avoid pbc errors..

ats = read(args.input, index=":")

cell_center = np.diag(ats[0].get_cell())/2.0
print cell_center

print 'Full Cell CrackPos: ', crack_pos
crack_pos = crack_pos - qm_com + cell_center
print 'Shift CrackPos', crack_pos

print 'QM Radius', qm_radius
print len(ats)
for at in ats[args.start:]:
#   mark quantum atoms
    x,y,z = at.positions.T
    radius1 = np.sqrt((x - crack_pos[0])**2 + (y-crack_pos[1])**2 + (z-crack_pos[2])**2)
    qm_region_mask = np.array(map(int, (radius1 <= qm_radius)))
    print sum(qm_region_mask)
    at.new_array("qm_atoms", qm_region_mask, dtype=float)
    at.set_array("qm_atoms", qm_region_mask)
    write_xyz("forces.xyz", at, append=True)


