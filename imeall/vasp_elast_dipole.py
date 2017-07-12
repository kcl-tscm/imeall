from quippy import Atoms, set_fortran_indexing 
from ase.calculators.vasp import Vasp
from imeall.calc_elast_dipole import ElasticDipole

set_fortran_indexing(False)

ed = ElasticDipole()

#We require the position vector of the defect atom in the cell.
ats_orig = Atoms('../defect_cell_relaxed.xyz')
h_list = [at for at in ats_orig if at.number==1]
defect = h_list[0]

vs = Vasp(restart=True)
ats = vs.get_atoms()
forces = vs.read_forces(ats)
ats.set_momenta(forces)
ats.write('test.xyz')

print ed.compute_vacancy_dipole(defect, ats, forces=forces)
