import re
import argparse
import numpy as np
from quippy import Atoms, set_fortran_indexing

set_fortran_indexing(False)

class Particle(object):
    def __init__(self, site, symbol, energy):
        self.symbol = symbol
        self.site = site
        self.energy = energy

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_file')
args = parser.parse_args()

h_line_re = re.compile(r'\[\s+([\-0-9\.]+)\s+([\-0-9\.]+)\s+([\-0-9\.]+)\s?\s?\]\s+([\-0-9\.]+)\s+([\-0-9\.]+)', re.S)

with open(args.input_file) as f:
    lines = h_line_re.findall(f.read())

interstitials = []
for line in lines:
    site = map(float, [line[0], line[1], line[2]])
    energy = float(line[4])
    interstitials.append(Particle(site,'H', energy))

interstitials.sort(key=lambda x:x.energy)
for int_ in interstitials:
    print int_.site, int_.energy

ats = Atoms('output.xyz')
num_fe = len(ats)
for int_ in interstitials:
    ats.add_atoms(np.array(int_.site), 1)

for i in range(len(ats)):
    ats.id[i] = i

ats.add_property('locen', 0.0)

for i, at in enumerate(interstitials):
    ats.properties['locen'][i+num_fe] = at.energy

ats.write('h_energetics.xyz')
