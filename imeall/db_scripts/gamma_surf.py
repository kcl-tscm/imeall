import os
import numpy as np
from   quippy            import set_fortran_indexing
import ase.units as units
from   ase.lattice.cubic import Diamond, BodyCenteredCubic
from   ase.optimize      import LBFGS
from   ase.constraints   import UnitCellFilter, StrainFilter, FixedLine, FixAtoms
from   ase.calculators.vasp import Vasp
from   ase.io.vasp       import write_vasp

from quippy import Atoms, Potential

set_fortran_indexing(False)

def gamma_surf(h_pos=np.array([1.41, 1.500, 22.48])):
    """
    :method:`gamma_surf` generates a set of directories with the upper
    half unit cell displaced along a certain distance
    along a particular lattice vector. Hydrogens can be added by
    setting vector in h_pos.

    TODO:
      h_pos should be a list of vectors and atom type for this to be more general.
    """

    vasp_exe = '/projects/SiO2_Fracture/iron/vasp.bgq'
    crack_geom = {'cleavage_plane'    : (1,1,0),
                  'crack_front'       : (1,-1,0),
                  'crack_direction'   : (0,0,1)}

    crack_direction = crack_geom['crack_direction']
    crack_front     = crack_geom['crack_front']
    cleavage_plane  = crack_geom['cleavage_plane']
    fe_unit = BodyCenteredCubic(directions=[crack_direction, crack_front, cleavage_plane],
                                size=(1,1,1), symbol='Fe', pbc=(1,1,1),
                                latticeconstant=2.83)
    nunits  = 8
    fe_bulk = BodyCenteredCubic(directions=[crack_direction, crack_front, cleavage_plane],
                                size=(2,2,nunits), symbol='Fe', pbc=(1,1,1),
                                latticeconstant=2.83)
    fe_unit = Atoms(fe_unit)
    fe_bulk = Atoms(fe_bulk)

    ycut    = 5.0 + float(nunits)/2.0*fe_unit.lattice[2,2]
    fe_bulk.center(vacuum=5.0, axis=2)
    print 'lattice', fe_bulk.lattice[1,1]
    print 'ycut', ycut

    a1 =  fe_unit.lattice[0,0]
    a2 =  fe_unit.lattice[1,1]

    POT_DIR = os.environ['POTDIR']
    eam_pot = os.path.join(POT_DIR, 'PotBH.xml')
    r_scale = 1.00894848312
    #eam_pot = os.path.join(POT_DIR, 'iron_mish.xml')
    #r_scale = 1.0129007626
    pot     = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(r_scale), param_filename=eam_pot)
    fe_bulk.set_calculator(pot)
    print 'Bulk Energy', fe_bulk.get_potential_energy()
    refen = fe_bulk.get_potential_energy()
    A = fe_bulk.get_volume()/fe_bulk.lattice[2,2]

    vasp_args=dict(xc='PBE', amix=0.01, amin=0.001, bmix=0.001, amix_mag=0.01, bmix_mag=0.001,
                   kpts=[8, 8, 1], kpar=32, lreal='auto', ibrion=2, nsw=40, nelmdl=-15, ispin=2,
                   nelm=100, algo='VeryFast', npar=8, lplane=False, lwave=False, lcharg=False, istart=0,
                   voskown=1, ismear=1, sigma=0.1, isym=2) # possibly try iwavpr=12, should be faster if it works

    dir_name  = 'b111shiftsym'
    #dir_name  = 'b001shiftsym'
    f         = open('_patdon123{0}H1.dat'.format(dir_name),'w')
    WRITEVASP = True
    a0        = 2.83
    for inum, ashift in enumerate(np.arange(0, 1.10, 0.10)):
        try:
            os.mkdir('{0}{1}'.format(dir_name, ashift))
        except:
            pass
        print 'Directory already exists'

        os.chdir('{0}{1}'.format(dir_name, ashift))
        fe_shift = fe_bulk.copy()
        fe_shift.set_calculator(pot)
        for at in fe_shift:
            if at.position[2] > ycut:
            #[00-1](110)
            #  at.position += ashift*np.array([-a1,0,0])
            #[1-11](110)
                at.position += 0.5*ashift*np.array([a1, a2,0])
    #  print >> fil1, ashift, (units.m**2/units.J)*(fe_shift.get_potential_energy()-refen)/A
        line=[]
        for at in fe_shift:
            line.append( FixedLine(at.index, (0,0,1)) )
    #Now add Hydrogen
    # fe_shift.add_atoms(np.array([0.53*a0, 0.53*a0, 21.0+a0/2]),1)
    # at relaxed position:
        fe_shift.add_atoms(h_pos,1)
        fix_atoms_mask = [at.number==1 for at in fe_shift]
        fixedatoms = FixAtoms(mask=fix_atoms_mask)
        fe_shift.set_constraint(line+[fixedatoms])
        opt   = LBFGS(fe_shift)
        opt.run(fmax=0.1, steps=1000)
        if inum==0:
            print 'Setting Reference Energy'
            refen = fe_shift.get_potential_energy()
        print >> f, ashift, (units.m**2/units.J)*(fe_shift.get_potential_energy()-refen)/A
        fe_shift.write('feb{0}.xyz'.format(ashift))
        if WRITEVASP:
            vasp = Vasp(**vasp_args)
            vasp.initialize(fe_shift)
            write_vasp('POSCAR', vasp.atoms_sorted,
                       symbol_count=vasp.symbol_count,
                       vasp5=True)
            vasp.write_incar(fe_shift)
            vasp.write_potcar()
            vasp.write_kpoints()
        os.chdir('../')
    f.close()

if __name__=="__main__":
    gamma_surf()
