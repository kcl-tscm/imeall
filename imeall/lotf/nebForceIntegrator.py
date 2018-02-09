# -*- coding: utf-8 -*-

"""Modified ase/neb.py to allow energy-free barrier climbing if desired"""

import threading
from math import sqrt

import numpy as np
from scipy.interpolate import spleval,interp1d
from scipy.integrate import cumtrapz
import ase.parallel as mpi
from ase.calculators.calculator import Calculator
from ase.io import read
from ase.optimize import BFGS,FIRE
from ase.utils.geometry import find_mic

class NEB:
    def __init__(self, images, k=0.1, force_only=False, modified_neb=True,
                 climb=False, parallel=False, world=None):
        """
        Nudged elastic band.
        images: list of Atoms objects
            Images defining path from initial to final state.
        k: float or list of floats
            Spring constant(s) in eV/Ang.  One number or one for each spring.
        climb: bool
            Use a climbing image (default is no climbing image).
        parallel: bool
            Distribute images over processors.
        force_only: bool
            Omit calls to get_potential_energy() function.
        climb: bool
            When ~force_only perform force integration to find climbing knot
        modified_neb = bool
            Use Henklemans kink-reduction NEB force
        """
        self.images = images
        self.climb = climb
        self.parallel = parallel
        self.natoms = len(images[0])
        self.nimages = len(images)
        self.emax = np.nan

        """ New members for spline integration"""
        self.force_only = force_only
        self.forces = np.zeros((self.nimages, self.natoms, 3))
        self.energies = np.zeros(self.nimages)
        self.modified_neb = modified_neb

        if isinstance(k, (float, int)):
            k = [k] * (self.nimages - 1)
        self.k = list(k)

        if world is None:
            world = mpi.world
        self.world = world


        if parallel:
            assert world.size == 1 or world.size % (self.nimages - 2) == 0

    def interpolate(self, method='linear', mic=True):
        interpolate(self.images, mic)
        if method == 'idpp':
            self.idpp_interpolate(traj=None, log=None, mic=mic)

    def idpp_interpolate(self, traj='idpp.traj', log='idpp.log', fmax=0.1,
                         optimizer=BFGS, mic=False):
        d1 = self.images[0].get_all_distances(mic=mic)
        d2 = self.images[-1].get_all_distances(mic=mic)
        d = (d2 - d1) / (self.nimages - 1)
        old = []
        for i, image in enumerate(self.images):
            old.append(image.calc)
            image.calc = IDPP(d1 + i * d, mic=mic)
        opt = optimizer(self, trajectory=traj, logfile=log)
        opt.run(fmax=fmax)
        for image, calc in zip(self.images, old):
            image.calc = calc

    def get_positions(self):
        positions = np.empty(((self.nimages - 2) * self.natoms, 3))
        n1 = 0
        for image in self.images[1:-1]:
            n2 = n1 + self.natoms
            positions[n1:n2] = image.get_positions()
            n1 = n2
        return positions

    def set_positions(self, positions):
        n1 = 0
        for image in self.images[1:-1]:
            n2 = n1 + self.natoms
            image.set_positions(positions[n1:n2])
            n1 = n2

            # Parallel NEB with Jacapo needs this:
            try:
                image.get_calculator().set_atoms(image)
            except AttributeError:
                pass

    def force_integrator(self,force=False,spline_points=100):
        tot_d = np.linalg.norm(find_mic(self.images[-1].get_positions() -
                            self.images[0].get_positions(),
                            self.images[0].get_cell(), self.images[0].pbc)[0])
        knot_dist = np.zeros(self.nimages)
        _U = np.zeros((self.nimages,3*self.natoms))
        _F = np.zeros((self.nimages,3*self.natoms))
        _U[0] = self.images[0].get_positions().reshape(3* self.natoms)
        _F[0] = self.forces[0].reshape(3*self.natoms)
        for _i,img in enumerate(self.images[1:]):
            _d = find_mic(img.get_positions() -
                        self.images[0].get_positions(),
                        self.images[0].get_cell(), self.images[0].pbc)[0]
            _d_b = np.linalg.norm(find_mic(img.get_positions() -
                        self.images[-1].get_positions(),
                        self.images[-1].get_cell(), self.images[-1].pbc)[0])
            knot_dist[_i+1] = .5*np.linalg.norm(_d)/tot_d + .5-.5*_d_b/tot_d
            _d += self.images[0].get_positions()
            _U[_i+1] = _d.reshape(3* self.natoms)
            _F[_i+1] = self.forces[_i+1].reshape(3*self.natoms)

        print 'knot_dist = ', list(knot_dist)
        U_spl = interp1d(knot_dist, _U.T, 'cubic')
        F_spl = interp1d(knot_dist, _F.T, 'cubic')

        #r = np.linspace(0., 1.+1./(spline_points-1), spline_points+1, endpoint=False)
        r = np.linspace(0., 1.0, spline_points, endpoint=True)
        dU_spl = U_spl._spline.derivative()
        dU = dU_spl(r)
        F = F_spl(r).T
        #dU = spleval(U_spl._spline, r, deriv=1)
        #F = spleval(F_spl._spline, r, deriv=0)

        splined_force = (dU * F).sum(axis=1)
        splined_energy = -cumtrapz(splined_force, r)

        if force:
            return r[:-1], splined_energy, splined_force[:-1]
        else:
            return r[:-1], splined_energy

    def get_forces(self):
        """Evaluate and return the forces."""
        images = self.images
        forces = np.empty(((self.nimages - 2), self.natoms, 3))
        energies = np.empty(self.nimages-2)
        if not self.parallel:
            for i in range(1, self.nimages - 1):
                if self.force_only:
                    energies[i - 1] = 0.
                else:
                    energies[i - 1] = images[i].get_potential_energy()
                forces[i - 1] = images[i].get_forces()

        elif self.world.size == 1:
            def run(image, energies, forces):
                if self.force_only:
                    energies[:] = 0.
                else:
                    energies[:] = image.get_potential_energy()
                forces[:] = image.get_forces()
            threads = [threading.Thread(target=run,
                                        args=(images[i],
                                              energies[i - 1:i],
                                              forces[i - 1:i]))
                       for i in range(1, self.nimages - 1)]
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join()
        else:
            # Parallelize over images:
            i = self.world.rank * (self.nimages - 2) // self.world.size + 1
            try:
                if self.force_only:
                    energies[i - 1] = 0.
                else:
                    energies[i - 1] = images[i].get_potential_energy()
                forces[i - 1] = images[i].get_forces()
            except:
                # Make sure other images also fail:
                error = self.world.sum(1.0)
                raise
            else:
                error = self.world.sum(0.0)
                if error:
                    raise RuntimeError('Parallel NEB failed!')

            for i in range(1, self.nimages - 1):
                root = (i - 1) * self.world.size // (self.nimages - 2)
                self.world.broadcast(energies[i - 1:i], root)
                self.world.broadcast(forces[i - 1], root)

        for i in range(1, self.nimages - 1):
            self.forces[i] = forces[i-1].copy()
            self.energies[i] = energies[i-1].copy()

        if self.force_only and self.climb:
            self.energies = self.get_potential_energies()
            energies = self.energies[1:-1]

        if self.force_only == True and self.climb == False:
            imax = self.nimages // 2  # just to avoid non-assigned error
        else:
            imax = 1 + np.argsort(energies)[-1]

        self.emax = energies[imax - 1]

        # Backwards tangent
        tangent1 = find_mic(images[1].get_positions() -
                                images[0].get_positions(),
                                images[0].get_cell(), images[0].pbc)[0]
        for i in range(1, self.nimages - 1):
            # Forwards tangent
            tangent2 = find_mic(images[i + 1].get_positions() -
                                images[i].get_positions(),
                                images[i].get_cell(),
                               images[i].pbc)[0]
            if i < imax:
                tangent = tangent2
            elif i > imax:
                tangent = tangent1
            else:
                tangent = tangent1 + tangent2
            if self.modified_neb:
                cos_t = np.vdot(tangent1, tangent2) / \
                    np.linalg.norm(tangent1) / np.linalg.norm(tangent2)
                f_c_t = .5 * (1. + np.cos(np.pi * cos_t))
            else:
                f_c_t = 0.
            tt = np.vdot(tangent, tangent)
            f = forces[i - 1]
            ft = np.vdot(f, tangent)
            if i == imax and self.climb:
                f -= 2 * ft / tt * tangent
            else:
                # -= F.T/|T| * T/|T| perp. MD force
                f -= ft / tt * tangent
                # f -= k(|T2| - |T1|) * (T/|T|) ~ ll spring force
                f -= (self.k[i-1]*np.linalg.norm(tangent1)-self.k[i] * \
                        np.linalg.norm(tangent2)) * tangent / np.sqrt(tt)
                if self.modified_neb:
                    f -= np.vdot(tangent1 * self.k[i - 1] - tangent2 * \
                                 self.k[i], tangent) / tt * tangent * (-f_c_t)
                    f -= (tangent1 * self.k[i - 1] - \
                          tangent2 * self.k[i]) * f_c_t

            tangent1 = tangent2.copy()
        return forces.reshape((-1, 3))

    def dump(self,file_name):
        al, spl_E, spl_F = self.force_integrator(force=True)
        al = al.reshape((-1,1))
        spl_E = spl_E.reshape((-1,1))
        spl_F = spl_F.reshape((-1,1))
        np.savetxt(file_name,np.hstack((al,spl_E,spl_F)))

    def get_potential_energies(self):
        _tmp_E = np.zeros(self.nimages)
        al, spl_E = self.force_integrator()
        _tmp_E[0]=0.
        _tmp_E[-1] = spl_E[-1]
        tot_d = np.linalg.norm(find_mic(self.images[-1].get_positions() -
                            self.images[0].get_positions(),
                            self.images[0].get_cell(), self.images[0].pbc)[0])
        for i in range(1,self.nimages - 1):
            _d = find_mic(self.images[i].get_positions() -
                        self.images[0].get_positions(),
                        self.images[0].get_cell(), self.images[0].pbc)[0]
            _d_b = np.linalg.norm(find_mic(self.images[i].get_positions() -
                        self.images[-1].get_positions(),
                        self.images[-1].get_cell(), self.images[-1].pbc)[0])
            _dist = .5*np.linalg.norm(_d)/tot_d + .5-.5*_d_b/tot_d
            _ind = min(int(len(spl_E)*_dist),len(spl_E)-1)
            _tmp_E[i] = spl_E[_ind]
        return _tmp_E

    def get_potential_energy(self, force_consistent=False):
        return self.emax

    def __len__(self):
        return (self.nimages - 2) * self.natoms

class IDPP(Calculator):
    """Image dependent pair potential.
    See:
        Improved initial guess for minimum energy path calculations.
        Søren Smidstrup, Andreas Pedersen, Kurt Stokbro and Hannes Jónsson
        Chem. Phys. 140, 214106 (2014)
    """

    implemented_properties = ['energy', 'forces']

    def __init__(self, target, mic):
        Calculator.__init__(self)
        self.target = target
        self.mic = mic

    def calculate(self, atoms, properties, system_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        P = atoms.get_positions()
        d = []
        D = []
        for p in P:
            Di = P - p
            if self.mic:
                Di, di = find_mic(Di, atoms.get_cell(), atoms.get_pbc())
            else:
                di = np.sqrt((Di**2).sum(1))
            d.append(di)
            D.append(Di)
        d = np.array(d)
        D = np.array(D)

        dd = d - self.target
        d.ravel()[::len(d) + 1] = 1  # avoid dividing by zero
        d4 = d**4
        e = 0.5 * (dd**2 / d4).sum()
        f = -2 * ((dd * (1 - 2 * dd / d) / d**5)[..., np.newaxis] * D).sum(0)
        self.results = {'energy': e, 'forces': f}


def interpolate(images, mic=False):
    """Given a list of images, linearly interpolate the positions of the
    interior images."""
    pos1 = images[0].get_positions()
    pos2 = images[-1].get_positions()
    d = pos2 - pos1
    if mic:
        d = find_mic(d, images[0].get_cell(), images[0].pbc)[0]
    d /= (len(images) - 1.0)
    for i in range(1, len(images) - 1):
        images[i].set_positions(pos1 + i * d)
        # Parallel NEB with Jacapo needs this:
        try:
            images[i].get_calculator().set_atoms(images[i])
        except AttributeError:
            pass
