import os
import sys
from models import GBMaintenance
import argparse

#"""""""""""""""""""""""""""""""""""""""#
#"Removes xyz files of different types."#
#"""""""""""""""""""""""""""""""""""""""#

gbm    = GBMaintenance()
parser = argparse.ArgumentParser()
parser.add_argument("-p", "--prefix", help = "Subsequent commands will act on all \
                     subdirectories with first characters matching prefix.", required=True)
parser.add_argument("-n", "--dryrun", help="If True print files which will be deleted with their size \
                                            If False then the files are truly deleted.", action ='store_true')
parser.add_argument("-rt", "--remove_type", help="rbt, structs, or traj, remove the base xyz and xyz.idx files for trajectories,\
                                                  base structure files or all subdirectories with rigibd body translation.")
args        = parser.parse_args()
prefix      = args.prefix
dryrun      = args.dryrun
remove_type = args.remove_type

for thing in os.listdir('./'):
    if os.path.isdir(thing) and thing[:len(prefix)] == prefix:
        gbm.remove_xyz(thing, dryrun=dryrun, remove_type=remove_type)
