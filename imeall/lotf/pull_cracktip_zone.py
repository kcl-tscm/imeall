import argparse
import numpy as np

from ase.io.xyz import write_xyz
from quippy import AtomsReader
from quippy import set_fortran_indexing

set_fortran_indexing(False)
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--radius", default=40.0, type=float, help="radius around cracktip to cutout (default 40\A)")
parser.add_argument("-i", "--input", default="crack_traj.xyz", help="trajectory file to cut crack tip from.")
parser.add_argument("-o", "--output", default="cracktip_zone.xyz", help="trajectory file to write cracktip region to.")

args = parser.parse_args()

def append(ats, rr, initial_crack, output_file='joined.xyz'):
    num_images = len(ats)
    for i, at in enumerate(ats):
        print i+1, '/', num_images, initial_crack
        fixed_mask = (np.sqrt(map(sum, map(np.square, at.positions[:,0:3]-initial_crack[0:3]))) <= rr)
        cl = at.select(fixed_mask)
        write_xyz(output_file, cl, append=True)

ats = AtomsReader(args.input)
initial_crackpos = np.array(ats[0].params['CrackPos'])
append(ats, args.radius, initial_crackpos, output_file=args.output)

