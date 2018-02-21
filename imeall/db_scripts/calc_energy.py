import os
import sys
from quippy import Potential, AtomsReader
import ase.units as units



class PotentialParameters(object):
    def __init__(self):
        self.name = 'Potential Parameters'

    def gs_ener_per_atom(self):
        eperat = {'Fe_Mendelev.xml' :-4.12243503431,
                  'Fe_Mendelev_Untruncated.xml' :-4.12243499408,
                  'Fe_Dudarev.xml'  : -4.31608690638,
                  'iron_mish.xml'   :-4.28000356875,
                  'Fe_Ackland.xml'  : -4.01298226805}
        return eperat

def calc_e_gb(at, E_bulk):
    cell = at.get_cell()
    A    = cell[0,0]*cell[1,1]
    E_gb = (at.get_potential_energy()-(at.n*(E_bulk)))/(2.*A)
    print A
    print at.get_potential_energy()
    print E_gb, 'eV/A^2'
    E_gb = 16.02*(at.get_potential_energy()-(at.n*(E_bulk)))/(2.*A)
    print E_gb, 'J/m^2'
    return E_gb

potparam = PotentialParameters()
ener_bulk = potparam.gs_ener_per_atom()

POT_DIR='/Users/lambert/pymodules/imeall/imeall/potentials'

at = AtomsReader(sys.argv[1])
at = at[-1]

pot_string = 'iron_mish.xml'
print ''
print '\t POTENTIAL: ', pot_string
print ''
eam_pot = os.path.join(POT_DIR, pot_string)
pot_1   = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(1.0129007626), param_filename=eam_pot)
at.set_calculator(pot_1)
calc_e_gb(at, ener_bulk[pot_string])

pot_string = 'Fe_Mendelev.xml'
print ''
print 'POTENTIAL: ', pot_string
print ''
eam_pot = os.path.join(POT_DIR, pot_string)
pot_2   = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(1.00894848312), param_filename=eam_pot)
at.set_calculator(pot_2)
calc_e_gb(at, ener_bulk[pot_string])

pot_string = 'Fe_Mendelev_Untruncated.xml'
print ''
print '\t POTENTIAL:', pot_string
print ''
eam_pot = os.path.join(POT_DIR, pot_string)
pot_2   = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(1.00894848312), param_filename=eam_pot)
at.set_calculator(pot_2)
calc_e_gb(at, ener_bulk[pot_string])

pot_string = 'Fe_Ackland.xml'
print ''
print '\t POTENTIAL: ', pot_string
print ''
eam_pot = os.path.join(POT_DIR, pot_string)
pot_2   = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(1.00894185389), param_filename=eam_pot)
at.set_calculator(pot_2)
calc_e_gb(at, ener_bulk[pot_string])

pot_string = 'Fe_Dudarev.xml'
print ''
print 'POTENTIAL:', pot_string
print ''
eam_pot = os.path.join(POT_DIR, pot_string)
pot_2   = Potential('IP EAM_ErcolAd do_rescale_r=T r_scale={0}'.format(1.00894185389), param_filename=eam_pot)
at.set_calculator(pot_2)
calc_e_gb(at, ener_bulk[pot_string])
