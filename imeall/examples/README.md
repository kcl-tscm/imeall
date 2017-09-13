checkout_structures:
  Demonstrates how to checkout all minimum energy structures for [001] alphaFe 
  grain boundaries for subsequent analysis.

query_test:
  Demonstrates a simple query to the database to pull the 
  minimum energy values for a user chosen potential 
  (or list of potentials), and orientation axis.
    $python query_struct.py -p 'PotBH.xml' -or 0,0,1

elastic_dipole_eam:
  Calculate the elastic dipole tensor for a hydrogen atom
  in bulk Fe.
    $python ~/imeall/imeall/calc_elast_dipole.py -ct EAM -i bcc_h.xyz -rc

dissolution_test:
  Calculate dissolution energy for hydrogen in bulk BCC Fe.

