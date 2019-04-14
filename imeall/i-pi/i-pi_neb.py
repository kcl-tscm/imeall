import os
import sys
import glob
import numpy as np

from ase import io
from ase.neb import fit0
from ase.io.extxyz import write_xyz
from ase import units

import matplotlib.pyplot as plt

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


def calc_neb_barriers(prefix, images):
    R = [ats.positions for ats in images]
    E = [ats.get_potential_energy() for ats in images]
    F = [atoms.get_forces() for atoms in images]
    A = images[0].cell
    pbc = images[0].pbc
    s, E, Sfit, Efit, lines = fit0(E, F, R, A, pbc)

#Normalizes path length to 1.
    s = np.array(s)
    norm = s.max()
    s /= norm
    Sfit /= norm

    fig = plt.figure()
    ax = fig.add_subplot(111)
    for x, y in lines:
        x /= norm
        ax.plot(x, y, '-C0')

    ax.plot(Sfit, Efit, 'C0-')
    Ef = max(Efit) - E[0]
    ax.plot(s, E, "oC0", label="%s: $E_f = %.4f$ (meV)"%("Interpolated Energy Barrier", Ef*1000))
    ax.legend(loc="best")
    ax.grid(True, linestyle="dashed")
    np.savetxt(prefix + "_s_E.txt", np.array([s, E]))
    np.savetxt(prefix + "_fit.txt", np.array([Sfit, Efit]))
    np.savetxt(prefix + "_lines.txt", np.vstack(lines))
    fig.savefig(prefix + "_barrier.eps")

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
        for step, row in enumerate((ener_array.T)[1:,:]):
            print >>f, step, ' '.join(map(str, row))
    return ener_array


if __name__=="__main__":
    print "prefix:", sys.argv[1]
    prefix = sys.argv[1]
    positions = glob.glob("{prefix}.pos_*xyz".format(prefix=prefix))
    forces = glob.glob("{prefix}.for_*xyz".format(prefix=prefix))
    
    #sort according to numerical index
    func = lambda x: int(x.split("_")[1].split(".")[0])
    positions.sort(key = func)
    forces.sort(key=func)
    
    potentials = read_neb_energies(prefix)
    #zero and convert to eV
    print np.shape(potentials)
    if np.shape(potentials)[0] > 1:
        potentials[1:,:] -= potentials[1:,:].min()
        potentials[1:,:] *= units.Hartree
    else:
        print potentials[0][1:] 
        potentials[0][1:] -= potentials[:].min()
        potentials[0][1:] *= units.Hartree

    print "Removing old neb path"
    if os.path.isfile("neb_path.xyz"):
        os.remove("neb_path.xyz")
    
    images = []
    print "shape of potential array", np.shape(potentials)
    #print "Final image potentials:", potentials[-1,1:]
    if np.shape(potentials)[0] > 1:
        final_pots = (potentials)[-1,1:]
    else:
        final_pots = potentials[0][1:]
    print final_pots, potentials
    for potential, pos_file, force_file in zip(final_pots, positions, forces):
        print pos_file, force_file, potential
        ats = io.read(pos_file,index="-1")
        alat = 14.35
        ats.cell = np.array([[alat,0.0,0.0],[0.0,alat,0.0],[0.0,0.0,alat]])
        ats.set_pbc(True) 
        forces = io.read(force_file,index="-1")
        #if conversion required:
        #forces *= units.Hartree/units.Bohr
        #bit of a hack but the positions of the force images are the forces:
        calc = Calculator(ats,potential,forces=forces.get_positions())
        ats.set_calculator(calc)
        ats.get_forces()
        images.append(ats)
        write_xyz("neb_path.xyz",ats, append=True)
    
    calc_neb_barriers(prefix, images)

