python ../../calc_elast_dipole.py -ct EAM -i bcc_h.xyz -rc
python ~/pymodules/imeall/imeall/calc_elast_dipole.py -i fe_bcc_h.xyz -ct EAM
qsub qm_edip.pbs
