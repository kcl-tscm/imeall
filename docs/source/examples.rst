Examples
========

Example scripts for interacting with the grain boundary database. Located in the ``imeall/examples/``
directory.

query_test:
-----------
Demonstrates a simple query to the database to pull the 
minimum energy values for a user chosen potential 
(or list of potentials), and orientation axis.

    ``python query_struct.py -p 'PotBH.xml' -or 0,0,1``

elastic_dipole_test:
--------------------
Calculate the elastic dipole tensor for a hydrogen atom in bulk Fe.

    ``python ~/imeall/imeall/calc_elast_dipole.py -ct EAM -i bcc_h.xyz -rc``

dissolution_test:
-----------------
Calculate dissolution energy for hydrogen in bulk BCC Fe.

    ``python dissolution_energy.py``

gen_grainboundaries:
--------------------
Use quaternion algebra to generate a series of canonical grain boundaries.

    ``python generate_gb_example.py``

boundary_interstitials:
-----------------------
Generate an array of potential interstitial sites using the voronoi analysis
tools from the modifed `pylada defects <https://github.com/pylada/pylada-defects>'_ module
for a particular boundary and obtain the dissolution energetics.

    ``python calc_eseg_example.py``

checkout_structures:
--------------------
Demonstrates how to checkout all minimum energy structures for [001] alphaFe tilt symmetric
grain boundaries for subsequent analysis. (Requires a local version of the database).

    ``python checkout_dirs.py``

