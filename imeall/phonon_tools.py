import os 
import sys
import numpy as np
import shutil
from distutils import dir_util

from ase.io import vasp
from ase import io
from os_supp.pushd import pushd
from argparse import ArgumentParser
import json

#displacement in angstrom
#index of atom in vasp list

def find_indices(center=np.array([0.0, 0.0, 0.0]), cutoff=6.0):
    """
    Find a cluster around a point
    Args:
        cutoff(float) :: radius of cluster.
        center(np.array) :: center of cluster.
    Returns:
        cluster(quippy.Atoms) ::
    """
    ats = io.read('./unpert/POSCAR')
    print 'Find cluster around', center
    cluster = []
    for index, at in enumerate(ats):
        vec = at.position - center
        if np.linalg.norm(vec) <= cutoff:
            print index, at.position
            cluster.append([index, at.position.tolist(), at.symbol])
        else:
            for latt_vec in ats.cell:
#test minus a latt vec
                vec = at.position - latt_vec - center
                if np.linalg.norm(vec) <= cutoff:
                    print index, at.position
                    cluster.append([index, at.position.tolist(), at.symbol])
#Test plus a latt vec
                vec = at.position + latt_vec - center
                if np.linalg.norm(vec) <= cutoff:
                    print index, at.position
                    cluster.append([index, at.position.tolist(), at.symbol])
    print len(cluster), 'atoms in the cluster'
    indices = [c for c in cluster]
    with open('indices.json', 'w') as f:
        json.dump({'indices':indices}, f)
    return cluster

def gen_cluster_dirs(center=np.array([0.0,0.0,0.0]), cutoff=6.0):
    cluster = find_indices(center=center, cutoff=cutoff)
    #cluster is a list: [index, position]
    for at_line in cluster:
        print at_line[0], at_line[1]
        directions = [(0, 'x'), (1,'y'), (2,'z')]
        for cart in directions:
            gen_forceconstants(at_line[2], at_line[0], cart, displacement=0.1, hour=5, np=48)

def gen_forceconstants(species, index, cart, displacement=0.1, hour=5, np=48, template_dir='templates'):
    """
    generate directory and POSCAR perturbation

    Args:
        species(str): Name of perturbed atom.
        
    """
    target_dir = "{spec}_{index}_{cart}_disp".format(spec=species, index=index, cart=cart[1])
    dir_util.copy_tree(template_dir, target_dir)
    with pushd(target_dir) as ctx0:
#Fix the pbs string
        with open('vasp.sh','r') as f:
          pbs_str = f.read()
        pbs_str = pbs_str.format(hour=hour, index=index, np=np, species=species)
        with open('vasp.sh','w') as f:
          print >> f, pbs_str

        ats = vasp.read_vasp('POSCAR')
        print "before pert", ats[index].position[:]
        ats[index].position[cart[0]] += displacement
        print "after pert", ats[index].position[:]
        vasp.write_vasp("POSCAR", ats)

def pull_forces(species, index, cart):
    """
    pull forces from OUTCAR and save to file.
    """
    target_dir = "{spec}_{index}_{cart}_disp".format(spec=species, index=index, cart=cart)
    with pushd(target_dir) as ctx0:
        print os.getcwd()
        ats = vasp.read_vasp_out()
        forces = ats.get_forces()
        np.save("col_vec", forces)
        ats.write("{}.xyz".format(target_dir))

if __name__=="__main__": 
    parser = ArgumentParser()
    parser.add_argument("-g", "--gen_forces", help="compute force constants for a displacement",  action="store_true")
    parser.add_argument("-p", "--pull_forces", help="pull the forces and store the 3N vector of force couplings.",  action="store_true")
    parser.add_argument("-s", "--species", default="Fe", type=str)
    parser.add_argument("-c", "--cart", default="x", type=str)
    parser.add_argument("-d", "--displacement", help="atomic displacement distance in angstrom", default=0.005, type=float)
    parser.add_argument("-i", "--index", help="index of atom in list", type=int)
    parser.add_argument("-t", "--time", help="hours to run.", default=2)
    parser.add_argument("-n", "--np", help="number of processors.", default=48)
    parser.add_argument("-a", "--array", help="generate array of jobs", action='store_true')
    parser.add_argument("-f", "--template_dir", help="folder containing POSCAR, INCAR templates", default='templates')

    args = parser.parse_args()

    if args.gen_forces:
        ats = io.read('./{tmp_folder}/POSCAR'.format(tmp_folder=args.template_dir))
        for at in ats:
            directions = [(0, 'x'), (1,'y'), (2,'z')]
            for cart in directions:
                gen_forceconstants(at.symbol, at.index, cart, displacement=args.displacement, hour=args.time, np=args.np, template_dir=args.template_dir)

    if args.pull_forces:
        pull_forces(args.species, args.index, args.cart)

    if args.array:
        gen_cluster_dirs(center=np.array([-0.07774,-11.4519, 14.0087]), cutoff=8.0)

