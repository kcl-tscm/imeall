import sys
from ase import io

mv_file = open("new_xbs.mv",'w')
images = io.read(sys.argv[1],index=":")
#start increment from one because the movie file
#is assumed to be from an xyz file where the first
#image is already in the .bs file.
for image in images[1:]:
    print >> mv_file, "frame"
    print >> mv_file, ' '.join(map(str,image.get_positions().flatten()))
    print >> mv_file, ""

mv_file.close()
