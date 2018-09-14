import os 
import json
import argparse
import numpy as np
import ase.units as units

from ase import io
from ase.io import Trajectory
from ase.io.xyz import write_xyz
from ase import Atoms as aseAtoms
from ase.constraints import FixAtoms
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution

from distutils import spawn

from imeall import app
from imeall.lotf.ForceMixerCarver import ForceMixingCarvingCalculator

from matscipy.socketcalc import VaspClient, SocketCalculator

from simulate_crack import update_qm_region_context, fix_edges, set_qmmm_pot, pass_print_context,\
                           check_if_cracked_context, pass_trajectory_context

from quippy import Atoms, set_fortran_indexing
from quippy.io import AtomsWriter, AtomsReader
from quippy.lotf import LOTFDynamics, update_hysteretic_qm_region
from quippy.crack import get_strain, get_energy_release_rate,\
                             ConstantStrainRate, find_crack_tip_stress_field
from quippy.clusters import HYBRID_NO_MARK, HYBRID_ACTIVE_MARK
from quippy.potential import Potential, ForceMixingPotential
from quippy.system import verbosity_set_minimum, verbosity_to_str

set_fortran_indexing(False)

sim_T       = units.kB*350.0 # Simulation temperature
nsteps      = 6000             # Total number of timesteps to run for
timestep    = 0.5*units.fs    # Timestep (NB: time base units are not fs!)
cutoff_skin = 2.0*units.Ang    # Amount by which potential cutoff is increased

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--restart", action="store_true", help="If present job is restarted from input_file, and not thermalized")
parser.add_argument("-i", "--input_file", help="file to restart from", default="surf_traj.xyz")
args = parser.parse_args()

mpirun = spawn.find_executable('mpirun')
vasp = '/home/mmm0007/vasp/vasp.5.4.1/bin/vasp_std'

vasp_args = dict(xc='PBE', amix=0.22, amin=0.02, bmix=0.8, amix_mag=1.10, bmix_mag=1.00,
                 kpts=[4, 1, 4], kpar=4, lreal='auto', nelmdl=-15, ispin=2, prec='Accurate', ediff=1.e-4,
                 nelm=100, algo='VeryFast', lplane=False, lwave=False, lcharg=False, istart=0, encut=420,
                 maxmix=35, voskown=0, ismear=1, sigma=0.1, isym=0)
procs = 96
# need to have procs % n_par == 0
n_par = 1
if procs <= 8:
    n_par = procs
else:
    for _npar in range(2, int(np.sqrt(1.*procs))):
        if procs % int(_npar) == 0:
            n_par = procs // int(_npar)

vasp_client = VaspClient(client_id=0, npj=procs, ppn=1,
                         exe=vasp, mpirun=mpirun, parmode='mpi',
                         ibrion=13, nsw=1000000,
                         npar=8, **vasp_args)

qm_pot = Potential(calculator=SocketCalculator(vasp_client))

if not args.restart:
    atoms = io.read('POSCAR')
else:
    atoms = AtomsReader(args.input_file)[-1]

atoms.set_calculator(qm_pot)

#thermalize atoms
if not args.restart:
    MaxwellBoltzmannDistribution(atoms, 2.0*sim_T)

dynamics = VelocityVerlet(atoms, timestep)

def print_context(ats=atoms, dyn=dynamics):
    print 'steps, T', dyn.nsteps, ats.get_kinetic_energy()/(1.5*units.kB*len(ats))

def write_slab(a=atoms):
    write_xyz('surf_traj.xyz', a, append=True)

dynamics.attach(print_context, interval=1)
dynamics.attach(write_slab, interval=1)
dynamics.run(nsteps)

