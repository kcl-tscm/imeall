This example requires spglib and tess to be installed. Both are available from the pip package manager.

Directory contains .xyz structure file for 0013687130 grain 
boundary min energy structure and associated json files. To generate
and calculate the interstitial sites for hydrogen dissolution run:

  ``python interstitial_example.py``

The relaxed positions of the H interstitals and the energetics of the site are
stored in the `h_site_ener_stretch_0.0.txt` file:

  ``python plot_h_energy.py -i ./h_site_ener_stretch_0.0.txt``

Open the resulting structure file in ovito, add the color coding modifier to visualize the energetics.
The property on the atoms object is called `locen`:
  ``ovito h_energetics.xyz``

Run clean.py to remove intermediate files:
  ``python clean.py``
