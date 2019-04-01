import os
import sys
import shutil
import argparse

from ase import io
from ase.units import Bohr
from os_supp.pushd import pushd

#rescale and put into units of alat
#site_str += ('ATOM=%-2s POS=%12.8f %12.8f %12.8f\n        '%((sym,) + tuple(pos/alat)))
def atoms_to_site(atoms):
    """
    Take an ase:Atoms: object and write it in site.dat format
    for lm-tbe code.
    """
    nbas = len(atoms)
    alat = 5.42351614
    lat_conv = 1.0/alat/Bohr

    header_line = "% site-data vn=3.0 fast io=62 nbas={nbas} alat={alat} plat= 5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0".format(nbas=nbas, alat=alat)
    comment_line = "#                        pos                                   vel                                    eula                   vshft  PL rlx"
    filler_line = "0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.000000  0  111\n"
    site_str = ''
    for symbol, pos in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
    #io=62 site strings
        pos_str = ' '.join(map(lambda x: '{:12.8f}'.format(x), lat_conv*pos))
        site_str += '{symbol:2s} {pos}  {fill}'.format(symbol=symbol, pos=pos_str, fill=filler_line)

    with open("site.dat",'w') as f:
        print >> f, header_line
        print >> f, comment_line
        print >> f, site_str

def pos_to_xyz(struct_file):
    with open(struct_file) as f:
        struct_str = f.read()
    struct_lines = struct_file.split('\n')[1:-1]
    ats.write("conv.xyz")


def gen_multiple_images():
  at_images = io.read("restart_images.xyz", index=":")
  print "{} images in file".format(len(at_images))
  for image_num, atoms in enumerate(at_images):
#Create separate directory for each image along path
    try:
        os.mkdir("image_{}".format(image_num))
    except OSError:
        print "dir {}  already exists".format(image_num)

    with pushd("image_{}".format(image_num)) as ctx1:
        shutil.copy("../ctrl.dat","./ctrl.dat")
        atoms_to_site(atoms)

if __name__=='__main__':

    parser = parser.ArgumentParser()
    parser.add_argument("-i","--input",help="input format suffix",default="xyz")
    parser.add_argument("-o","--output",help="output format suffix",default="site")
    parser.add_argument("-s","--struct_file",help="input format suffix",required=True)

    if args.input=="xyz" and args.output=="site":
        ats = io.read(args.struct_file)
        atoms_to_site(ats)
    if args.input=="pos" and args.output=="xyz":
        pos_to_xyz(args.struct_file)
    else:
        print "No conversion method available"


