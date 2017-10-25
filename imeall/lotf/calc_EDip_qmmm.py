from __future__ import print_function
import os
import sys

import json 
import numpy as np

from ase.optimize import FIRE
from ase.optimize.precon import PreconFIRE, Exp
import argparse

from distutils import spawn

from imeall import app
from ForceMixerCarver import ForceMixingCarvingCalculator
from matscipy.socketcalc import VaspClient, SocketCalculator

from quippy import AtomsReader, AtomsWriter
from quippy import Potential, set_fortran_indexing

from imeall.calc_elast_dipole import find_h_atom

set_fortran_indexing(True)

class NEBAnalysis(object):
  def __init__(self):
    pass

  def save_barriers(self, images, neb, prefix=""):
    """ Saves the NEB results and plots energy barrier.
    Arguments:
      images :: 
      prefix ::
      neb :: :class:ase:`NEB`
    """
    if not os.path.exists('NEB_images'):
        os.mkdir('NEB_images')

    for i, image in enumerate(images):
        write("NEB_images/" + prefix + "ini_NEB_image%03i.xyz" % i, image)

    plot_dir = "NEB_plots"
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    R = [atoms.positions for atoms in images]
    E = neb.get_potential_energies()
    F = [atoms.get_forces() for atoms in images]
    A = images[0].cell
    pbc = images[0].pbc
    s, E, Sfit, Efit, lines = fit0(E, F, R, A, pbc)

    s = np.array(s)
    norm = s.max()
    s /= norm
    Sfit /= norm

    for x, y in lines:
        x /= norm
        ax.plot(x, y, '-C0')

    ax.plot(Sfit, Efit, 'C0-')
    Ef = max(Efit) - E[0]
    ax.plot(s, E, "oC0", label="%s: $E_f = %.4f$"%("QM/MM EAM result", Ef))
    ax.legend(loc="best")
    ax.grid(True, linestyle="dashed")

    np.savetxt(plot_dir + '/' + prefix + "_s_E.txt", np.array([s, E]))
    np.savetxt(plot_dir + '/' + prefix + "_fit.txt", np.array([Sfit, Efit]))
    np.savetxt(plot_dir + '/' + prefix + "_lines.txt", np.vstack(lines))

    fig.savefig(plot_dir + '/' + prefix + "_barrier.eps")
    return None


def build_neb_configurations():
  at_ini =  AtomsReader('gb_traj.xyz')[-1]

mpirun = spawn.find_executable('mpirun')
vasp = '/home/mmm0007/vasp/vasp.5.4.1/bin/vasp_std'

parser = argparse.ArgumentParser()
parser.add_argument("--buff", "-b", type=float, default=6.0)
parser.add_argument("--qm_radius", "-q", type=float, default=5.0)
parser.add_argument("--use_socket", "-u", type=bool, default=True)
parser.add_argument("--input", "-inp", help=".xyz input filename.")
args=parser.parse_args()

buff = args.buff
qm_radius = args.qm_radius

POT_DIR = os.path.join(app.root_path, 'potentials')
eam_pot = os.path.join(POT_DIR, 'PotBH.xml')
r_scale = 1.00894848312
mm_pot = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)

#disloc_fin, disloc_ini, bulk = sd.make_barrier_configurations(calculator=lammps, cylinder_r=cylinder_r)
gb_cell = AtomsReader('bcc_h.xyz')[-1]
defect = find_h_atom(gb_cell)
h_pos = defect.position

x, y, z = gb_cell.positions.T
radius1 = np.sqrt((x - h_pos[0])**2 + (y-h_pos[1])**2 + (z-h_pos[2])**2)

qm_region_mask = (radius1 < qm_radius)
qm_buffer_mask = (radius1 < qm_radius + buff)

print ("\nNumber of atoms in qm region of %.1f" % qm_radius +
                                "A : %i" % np.count_nonzero(qm_region_mask))

print ("together with the buffer of %.1f" % (qm_radius + buff ) +
                                "A %i" % np.count_nonzero(qm_buffer_mask))

magmoms=[2.6 for _ in range(np.count_nonzero(qm_buffer_mask))]
vasp_args = dict(xc='PBE', amix=0.01, amin=0.001, bmix=0.001, amix_mag=0.01, bmix_mag=0.001,
                 kpts=[1, 1, 1], kpar=1, lreal='auto', nelmdl=-15, ispin=2, prec='Accurate',
                 nelm=100, algo='VeryFast', lplane=False, lwave=False, lcharg=False, istart=0, 
                 magmom=magmoms, maxmix=30, #https://www.vasp.at/vasp-workshop/slides/handsonIV.pdf #for badly behaved clusters.
                 voskown=0, ismear=1, sigma=0.1, isym=2) # possibly try iwavpr=12, should be faster if it works

# parallel config.
procs = 24
kpts = [1, 1, 1]
# need to have procs % n_par == 0
n_par = 1
if procs <= 8:
  n_par = procs
else:
  for _npar in range(2, int(np.sqrt(1.*procs))):
    if procs % int(_npar) == 0:
      n_par = procs // int(_npar)


if args.use_socket:
  vasp_client = VaspClient(client_id=0, npj=procs, ppn=1,
                           exe=vasp, mpirun=mpirun, parmode='mpi',
                           ibrion=13, nsw=1000000,
                           npar=n_par, **vasp_args)
  qm_pot = SocketCalculator(vasp_client)
else:
  pass

qmmm_pot = ForceMixingCarvingCalculator(gb_cell, qm_region_mask,
                                        mm_pot, #mm_pot_mod, #for testing
                                        qm_pot,
                                        buffer_width=buff,
                                        pbc_type=[False, False, False])


def pass_trajectory_context(trajectory, dynamics):
  def traj_writer(dynamics):
      trajectory.write(dynamics.atoms)
  return traj_writer

#Relax atoms around the defect.
trajectory = AtomsWriter('gb_traj_ini.xyz')
opt = FIRE(gb_cell)
opt.attach(pass_trajectory_context(trajectory, opt), 1, opt)
gb_cell.set_calculator(qmmm_pot)
opt.run(fmax=1.0e-3)
gb_cell.write('dft_edip_H_relaxed.xyz')

#One shot calculation of forces with defect removed.
gb_cell = AtomsReader('dft_edip_H_relaxed.xyz')[-1]
#Remove defect
defect = find_h_atom(gb_cell)
gb_cell.remove_atoms([defect.index+1])
#Get forces
gb_cell.set_calculator(qmmm_pot)
gb_cell.get_forces()
gb_cell.write('dft_edip_forces.xyz')

sock_calc.shutdown()
