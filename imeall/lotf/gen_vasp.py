import os
import sys
import ase.units as units
from ase.calculators.vasp import Vasp
from quippy import bcc
#from atomsserver import QUIPClient, VaspClient
#SCRIPT to set the magnetic moments automatically in INCAR file.

x = bcc(2.83)
x.set_atoms(26)
x = x*(4,4,1)
mom = [3.0 for at in x]
print len(mom), len(x)
x.set_initial_magnetic_moments(mom)
print x.get_initial_magnetic_moments()
# VASP arguments
vasp_args=dict(xc='PBE', amix=0.01, amin=0.01, bmix=0.001, amix_mag=0.01, bmix_mag=0.0001,
               kpts=[1, 1, 6], kpar=4, lreal='auto', ibrion=13, nsw=1000000, nelmdl=-15, ispin=2,
               nelm=100, algo='VeryFast', npar=32, lplane=False, lwave=False, lcharg=False, istart=0,
               voskown=1, ismear=1, sigma=0.1, isym=0, magmom=x.get_initial_magnetic_moments()) # possibly try iwavpr=12, should be faster if it works
vc = Vasp(**vasp_args)
vc.initialize(x)
vc.write_incar(x)
