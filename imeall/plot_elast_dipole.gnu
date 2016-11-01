set xlabel 'Mean H Spacing (A)'
set ylabel 'Energy (J/m^2)'
pl './gb/H_en.dat'             u 2:6 w lp pt 6 lc 2 t 'GB',\
   './gb_strained/H_en.dat'    u 2:6 w lp pt 7 lc 2 t 'GB Strained',\
   './gb_compress/H_en.dat'    u 2:6 w lp pt 8 lc 2 t 'GB Compressed',\
   './bulk/H_en.dat'           u 2:6 w lp pt 6 lc -1 t 'Bulk',\
   './bulk_strained/H_en.dat'  u 2:6 w lp pt 7 lc -1 t 'Bulk Strained',\
   './bulk_compress/H_en.dat'  u 2:6 w lp pt 8 lc -1 t 'Bulk Compressed'
