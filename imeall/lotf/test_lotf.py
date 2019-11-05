import os
import sys
import json
import argparse
import numpy as np
import ase.units as units

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
from quippy.clusters import HYBRID_NO_MARK, HYBRID_ACTIVE_MARK, create_pos_or_list_centred_hybrid_region
from quippy.potential import Potential, ForceMixingPotential
from quippy.system import verbosity_set_minimum, verbosity_to_str

set_fortran_indexing(False)
VERBOSITY = 2
verbosity_set_minimum(VERBOSITY)
print verbosity_to_str(VERBOSITY)

extrapolate_steps = 4        # Number of steps for predictor-corrector
                             # interpolation and extrapolation

with open('crack_info.json','r') as f:
  crack_dict = json.load(f)

sim_T       = crack_dict['sim_T'] # Simulation temperature
nsteps      = 6000             # Total number of timesteps to run for
timestep    = 1.0*units.fs    # Timestep (NB: time base units are not fs!)
cutoff_skin = 2.0*units.Ang    # Amount by which potential cutoff is increased
                               # for neighbour calculations
tip_move_tol = 12.0            # Distance tip has to move before crack
                               # is taken to be running
strain_rate    = 0.0
traj_interval  = 80             # Number of time steps between interpolations
print_interval = 10             # time steps between trajectory prints 10 fs
param_file     = 'params.xml'   # Filename of XML file containing
                                # potential parameters
mm_init_args = 'IP EAM_ErcolAd do_rescale_r=T r_scale=1.008948312' # Classical potential
qm_init_args = 'TB DFTB'        # Initialisation arguments for QM potential

input_file   = 'crack.xyz'      # crack_slab
traj_file    = 'crack_traj.xyz' # Trajectory output file in (NetCDF format)

qm_inner_radius = 6.0
qm_outer_radius = 8.5


try:
  pot_dir      = os.environ['POTDIR']
except:
  sys.exit("PLEASE SET export POTDIR='path/to/potfiles/')")

pot_file     = os.path.join(pot_dir, 'PotBH.xml')
# Potential information
POT_DIR = os.path.join(app.root_path, 'potentials')
eam_pot = os.path.join(POT_DIR, 'PotBH.xml')
mpirun = spawn.find_executable('mpirun')
vasp = '/home/mmm0007/vasp/vasp.5.4.1/bin/vasp_std'

atoms = AtomsReader('traj_1p5ps.xyz')[-1]
strain_atoms = fix_edges(atoms)

crack_pos = atoms.info['CrackPos']
r_scale = 1.00894848312
mm_pot = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot, cutoff_skin=2.0)


qm_classical= False

if qm_classical:
    r_scale = 0.99894848312
    qm_pot = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot, cutoff_skin=2.0)
else:
    vasp_args = dict(xc='PBE', amix=0.01, amin=0.001, bmix=0.001, amix_mag=0.01, bmix_mag=0.001,
                     kpts=[1, 1, 8], kpar=6, lreal='auto', nelmdl=-15, ispin=2, prec='Accurate', ediff=1.e-4,
                     nelm=100, algo='VeryFast', lplane=False, lwave=False, lcharg=False, istart=0, encut=420,
                     maxmix=30, voskown=0, ismear=1, sigma=0.1, isym=0)
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
                             npar=4, **vasp_args)

    qm_pot = Potential(calculator=SocketCalculator(vasp_client))

#construct_buffer_use_only_heavy_atoms might be necessary to stop confusion with H!
qmmm_pot = ForceMixingPotential(pot1=mm_pot, pot2=qm_pot, atoms=atoms,
                                qm_args_str='carve_cluster=T cluster_periodic_z=T cluster_calc_connect=T '+ 
                                            'single_cluster=T cluster_vacuum=5.0 terminate=F print_clusters=F',
                                fit_hops=4,
                                carve_cluster=True,
                                lotf_spring_hops=3,
                                terminate=False,
                                calc_connect=True,
                                buffer_hops=3,
                                hysteretic_buffer=True,
                                cluster_vacuum = 5.0,
                                cluster_hopping = True,
                                single_cluster=True,
                                hysteretic_buffer_inner_radius=qm_inner_radius,
                                hysteretic_buffer_outer_radius=qm_outer_radius,
                                cluster_hopping_nneighb_only=True,
                                min_images_only=True)

atoms.cutoff = 3.0
atoms.calc_connect()
qm_list = update_hysteretic_qm_region(atoms, [], crack_pos,
                                      qm_inner_radius, qm_outer_radius)
qmmm_pot.set_qm_atoms(qm_list, atoms)
atoms.set_calculator(qmmm_pot)
del_keys = ['force','forces', 'forces0']
array_keys = atoms.arrays.keys()
prop_keys = atoms.properties.keys()
for key in del_keys:
    if key in array_keys:
        del(atoms.arrays[key])

#under construction
if False:
    np.random.seed(42)
    MaxwellBoltzmannDistribution(atoms, 2.0*sim_T)

check_force_error=False
dynamics = LOTFDynamics(atoms, timestep, extrapolate_steps, check_force_error=check_force_error)

def log_pred_corr_errors(dynamics, logfile):
    logline = '%s err %10.1f%12.6f%12.6f\n' % (dynamics.state_label,
                                           dynamics.get_time()/units.fs,
                                           dynamics.rms_force_error,
                                           dynamics.max_force_error)
    print logline
    logfile.write(logline)

if check_force_error:
    pred_corr_logfile = open('pred-corr-error.txt','w')
    dynamics.attach(log_pred_corr_errors, 1, dynamics, pred_corr_logfile)

# array to store time averaged stress field
avg_sigma = np.zeros((len(atoms), 3, 3))

def update_qm_region(atoms=dynamics.atoms):
    #crack_pos = find_crack_tip_stress_field(atoms, calc=mm_pot, avg_sigma=avg_sigma)
    #print crack_pos
    crack_pos = atoms.info['CrackPos']
    qm_list   = qmmm_pot.get_qm_atoms(atoms)
    qm_list   = update_hysteretic_qm_region(atoms, qm_list, crack_pos,
                                            qm_inner_radius, qm_outer_radius)
    qmmm_pot.set_qm_atoms(qm_list, atoms)

print "Initialising Dynamics"
dynamics.set_qm_update_func(update_qm_region)

def write_slab(dynamics=dynamics):
    if dynamics.state == LOTFDynamics.Interpolation:
        dynamics.atoms.set_array('avg_sigma', avg_sigma.reshape((len(atoms), 9)))
        write_xyz('crack_traj_lotf.xyz', dynamics.atoms, append=True)

def write_cluster(ats=atoms, qmmm_pot=qmmm_pot):
    qm_list = list((qmmm_pot.atoms.hybrid > 0).nonzero()[0])
    qm_ats = ats[qm_list]
    write_xyz('cluster.xyz', qm_ats, append=True)

dynamics.attach(write_slab, interval=1)
print 'Running Dynamics'
dynamics.run(nsteps)

