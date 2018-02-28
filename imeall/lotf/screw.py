import os
import numpy as np
from ase.lattice.cubic import BodyCenteredCubic
from ase.optimize      import BFGS, FIRE, LBFGS, MDMin, QuasiNewton

from ase import Atom
from ase import Atoms as aseAtoms

from quippy.io import AtomsWriter, AtomsReader, write
from quippy.structures import disloc_noam
from quippy import Potential, Atoms, write, supercell, calc_nye_tensor
from quippy import set_fortran_indexing
from quippy.farray import fzeros, frange, unravel_index, farray

#Create unit cell with the orientation:

class Dislocation(object):
    def __init__(self, x=[1,1,-2], y=[-1,1,0],z=[1,1,1], name='defect.xyz'):
        """
        Args:
        x,y,z: lattice vector directions.
        """
        self.x = x
        self.y = y
        self.z = z
        self.l_x =  200.
        self.l_y =  200.
        self.l_z =  200.
        self.name = name

    def gen_screw_dislocation(self):
        screw_slab_unit = BodyCenteredCubic(directions = [self.x, self.y, self.z],
                                            size = (1,1,1), symbol='Fe', pbc=(1,1,1),
                                            latticeconstant = 2.83)
        screw_slab_unit = Atoms(screw_slab_unit)
# set the supercell:
        n_x = int(self.l_x/screw_slab_unit.get_cell()[0,0])
        n_y = int(self.l_y/screw_slab_unit.get_cell()[1,1])
        n_z = 2
        screw_slab_super = supercell(screw_slab_unit, n_x, n_y, n_z)
        ref_slab = screw_slab_unit*(n_x,n_y,n_z)
        ref_slab.write('ref_slab.xyz')
#Burgers vector modulus:
        b = screw_slab_unit.get_cell()[2,2]
        brg_vec = b*np.array([0,0,1])
        print b

        vacuum = 70.
        screw_slab_super.lattice[0,0] += vacuum
        screw_slab_super.lattice[1,1] += vacuum

        disloc_l   = np.array([0,0,1])
        screw_slab_super.set_lattice(screw_slab_super.get_cell(), scale_positions=False)
        super_slab = screw_slab_super.copy()

        L_x  = screw_slab_super.lattice[0,0]
        L_y  = screw_slab_super.lattice[1,1]
        L_z  = screw_slab_super.lattice[2,2]
        core =  np.array([(L_x-vacuum)/2., (L_y-vacuum)/2., L_z/2.])

        screw_slab_super.set_cutoff(3.0)
        screw_slab_super.calc_connect()
        disloc_noam(screw_slab_super, core, disloc_l, brg_vec)
        screw_slab_super.info['core'] = core
        screw_slab_super.write('s{0}.xyz'.format(self.name))


    def gen_quad_smart(self, dis_type="easy"):
        alat = 2.82893
        a1 = (1./3.)*alat*np.array([-1,-1,2])
        a2 = (1.0/2.0)*alat*np.array([1,-1,0])
        a3 = alat*0.5*np.array([1, 1, 1])
        M = np.array([a1,a2,a3])
        print M
        latt_1 = alat*np.array([-1.0,-1.0,2.0])
        latt_2 = alat*np.array([1.0, 0.0,-1.0])
        basis_1 = alat*np.array([0.0,0.0,0.0])
        basis_2 = alat*np.array([0.0,0.0,1.0])
        basis_3 = alat*np.array([0.5, 0.5, 0.5])
        #basis_1 = alat*np.array([0.0,0.0,0.0])
        #basis_2 = alat*np.array([0.81625, 0.0, 0.57734])
        #basis_3 = alat*np.array([1.6329, 0.0, 0.288675])
        #latt_1 = np.array([6.93206,0.0,0.0])
        #latt_2 = np.array([-3.46603,2.00111,0.0])
        #basis_2 = alat*np.sqrt(3.)*0.5*np.array([1.0,1.0,1.0])
        print 'new latt 1', np.dot(np.linalg.inv(M.T), basis_1)
        print 'new latt 2', np.dot(np.linalg.inv(M.T), basis_2)


        #as integeres
        if dis_type == "easy":
            n = 15.
            m = 9.
            c1 = n*a1 - (1.0/(3.0*m))*a3
            c2 = (n/2.)*a1 + m*a2 + (0.5 - (1.0/6*m))*a3
            c3 = a3
        elif dis_type == "hard":
            n = 21.
            m = 13.
            c1 = n*a1 + (1.0/(3*m))*a3
            c2 = (n/2.)*a1 + m*a2 + (0.5 + 1./(6.*m))*a3
            c3 = a3
        #screw_slab_unit = BodyCenteredCubic(directions = [c1,c2,c3],
        #                                    size = (1,1,1), symbol='Fe', pbc=(1,1,1),
        #                                    latticeconstant = 2.83)
        screw_slab_unit = aseAtoms()
        for x in range(-30,30):
            for y in range(-30,30):
                latt = x*(latt_1) + y*(latt_2)
                Fe_1 = latt
                Fe_2 = latt + basis_2
                Fe_3 = latt + basis_3
                #latt_norm = np.linalg.norm(latt)
                #c1_norm = np.linalg.norm(c1)
                #c2_norm = np.linalg.norm(c2)
                #c1_arg = np.arccos((np.dot(c1,latt)/(c1_norm*latt_norm))) 
                #c2_arg = np.dot(c2,latt)/(np.linalg.norm(c2)*np.linalg.norm(latt))
                #c1_comp = np.linalg.norm(np.dot(c1,latt)*(c1_norm*latt_norm))
                #c2_comp = np.linalg.norm(np.dot(c2,latt)*(c2_norm*latt_norm))
                screw_slab_unit.append(Atom(symbol='Fe', position=Fe_1))
                screw_slab_unit.append(Atom(symbol='Fe', position=Fe_2))
                screw_slab_unit.append(Atom(symbol='Fe', position=Fe_3))

        screw_slab_unit.set_cell([c1,c2,c3])
        screw_slab_unit.write('screw.xyz')
        
    def gen_screw_dislocation_quad(self):
        # Easy core
        screw_slab_unit = BodyCenteredCubic(directions = [self.x, self.y, self.z],
                                            size = (1,1,1), symbol='Fe', pbc=(1,1,1),
                                            latticeconstant = alat)
        screw_slab_unit = Atoms(screw_slab_unit)
# set the supercell:
        n_x = int(18.0/screw_slab_unit.get_cell()[0,0])
        n_y = int(18.0/screw_slab_unit.get_cell()[1,1])
        n_z = 2

        ref_slab = screw_slab_unit*(n_x,n_y,n_z)
        ref_slab.write('ref_slab.xyz')
        screw_slab_super = supercell(screw_slab_unit, n_x, n_y, n_z)
        b = screw_slab_unit.get_cell()[2,2]
        brg_vec = b*np.array([0,0,1])
        disloc_l = np.array([0,0,1])
        super_slab = screw_slab_super.copy()
        L_x  = screw_slab_super.lattice[0,0]
        L_y  = screw_slab_super.lattice[1,1]
        L_z  = screw_slab_super.lattice[2,2]
        core = [(np.array([(L_x)/2., (L_y)/2., L_z/2.]), brg_vec),
                 np.array([(L_x)/2., (L_y)/2., L_z/2.]),
                 np.array([(L_x)/2., (L_y)/2., L_z/2.]),
                 np.array([(L_x)/2., (L_y)/2., L_z/2.])]
        

    def gen_edge_dislocation(self):
        screw_slab_unit = BodyCenteredCubic(directions = [self.x, self.y, self.z],
                                            size = (1,1,1), symbol='Fe', pbc=(1,1,1),
                                            latticeconstant = 2.83)
        screw_slab_unit = Atoms(screw_slab_unit)

        n_x = int(self.l_x/screw_slab_unit.get_cell()[0,0])
        n_y = int(self.l_y/screw_slab_unit.get_cell()[1,1])
        n_z = 1

        screw_slab_super = supercell(screw_slab_unit, n_x, n_y, n_z)
        screw_slab_unit.write('ref_slab.xyz')

        b = screw_slab_unit.get_cell()[0,0]
        brg_vec = b*np.array([1,0,0])

        vacuum = 70.
        screw_slab_super.lattice[1,1] += vacuum

        disloc_l = np.array([0,0,1])
        screw_slab_super.set_lattice(screw_slab_super.get_cell(), scale_positions=False)
        super_slab = screw_slab_super.copy()

        L_x  = screw_slab_super.lattice[0,0]
        L_y  = screw_slab_super.lattice[1,1]
        L_z  = screw_slab_super.lattice[2,2]
#Maintain Periodicity along x and z for an edge dislocation:
        core =  np.array([(L_x)/2., (L_y-vacuum)/2., (L_z)/2.])

        screw_slab_super.set_cutoff(3.0)
        screw_slab_super.calc_connect()
        disloc_noam(screw_slab_super, core, disloc_l, brg_vec)
        screw_slab_super.info['core'] = core
        screw_slab_super.write('e{0}.xyz'.format(self.name))

if __name__=='__main__':
    disloc = Dislocation()
    disloc.gen_quad_smart()
