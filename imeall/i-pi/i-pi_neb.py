import os
import sys
import glob
import numpy as np
from ase import io
from ase.io.extxyz import write_xyz

class Calculator(object):
    def __init__(self,atoms,energy=0.0,forces=None):
        self.results = {'energy': energy,
                'forces': forces, #np.zeros((len(atoms), 3)),
                'stress': np.zeros(6),
                'dipole': np.zeros(3),
                'charges': np.zeros(len(atoms)),
                'magmom': 0.0,
                'magmoms': np.zeros(len(atoms))}

    def get_potential_energy(self,atoms):
        return self.results['energy']

    def get_forces(self,atoms):
        return self.results['forces']

    def get_stress(self, atoms):
        return np.zeros(6)

    def get_magmoms(self,atoms):
        """Return the Magnetic moments."""
        return self.results['magmoms']
    
    def calculation_required(self, atoms):
        return False


def read_neb_energies(prefix):
    with open('{}.neb'.format(prefix),'r') as f:
        neb_str = f.read()
    ener_list = []
    for line in neb_str.split('\n')[:-1]:
        if line[0]=='#':
            continue
        else:
            ener_list.append(np.array(map(float, line.split())))
    ener_array = np.array(ener_list)
    with open('neb_energies.dat','w') as f:
        for step, row in enumerate(ener_array.T[1:]):
            print >>f, step, ' '.join(map(str, row))

prefix = sys.argv[1]
positions = glob.glob("{prefix}.pos_*xyz".format(prefix=prefix))
forces = glob.glob("{prefix}.for_*xyz".format(prefix=prefix))

if os.path.isfile("neb_path.xyz"):
    os.remove("neb_path.xyz")

#sort according to numerical index
func = lambda x: int(x.split("_")[1].split(".")[0])
positions.sort(key = func)
forces.sort(key=func)

read_neb_energies(prefix)

for pos_file, force_file in zip(positions, forces):
    print pos_file, force_file
    ats = io.read(pos_file,index="-1")
    forces = io.read(force_file,index="-1")
#bit of a hack but the positions of the force images are the forces:
    calc = Calculator(ats,energy=0.0,forces=forces.get_positions())
    ats.set_calculator(calc)
    ats.get_forces()
    write_xyz("neb_path.xyz",ats, append=True)

