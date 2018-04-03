import os
import sys
import glob
import json
import shutil

import numpy as np
from os_supp.pushd import pushd

def canonical_name_from_gb_json(gb_json):
    """
    Generate a canonical name from the json file parameters.
    Args:
        gb_json: top level json dictionary.
    """
    degrees = gb_json['angle']*(180.0/np.pi)
    #degrees = str(np.round(degrees,2)).replace('.','')
    degrees = '{:0.2f}'.format(degrees).replace('.','')
    if len(degrees) == 5:
        pass
    elif len(degrees) == 4:
        degrees = '0'+degrees
    elif len(degrees) == 3:
        degrees = '00'+degrees
    else:
        sys.exit('Bad angle string')
    or_axis = '_'.join(map(str, map(int, gb_json['orientation_axis'])))
    or_axis = or_axis.replace('-','m')
    bp = '_'.join(map(str, map(int, gb_json['boundary_plane'])))
    bp = bp.replace('-','m')
    canonical_name = or_axis +'_{}_'.format(degrees) + bp
    return canonical_name

def rename_canonical_dirs():
    """
    Rename canonical dirs just below top level directory.
    """
    canonical_dirs = glob.glob('{}*'.format(or_x))
    canonical_dirs = filter(os.path.isdir, canonical_dirs)
    for canonical_dir in canonical_dirs:  
    #update top level canonical dir name to new convention
        with open(canonical_dir+'/gb.json') as f:
            gb_json  = json.load(f)
            old_canonical_name = gb_json['gbid']
        new_canonical_name =  canonical_name_from_gb_json(gb_json)
        print canonical_name_from_gb_json(gb_json)
        shutil.move(canonical_dir, new_canonical_name)
        with pushd(new_canonical_name) as ctx1:
            shutil.move(old_canonical_name+'.png', new_canonical_name+'.png')
            shutil.move(old_canonical_name+'.xyz', new_canonical_name+'.xyz')
            shutil.move('csl_'+old_canonical_name+'.svg', 'csl_'+new_canonical_name+'.svg')
            with open('gb.json') as f:
                gb_json = json.load(f)
                old_gbid = gb_json['gbid']
                gb_json['old_gbid'] = old_gbid
                gb_json['gbid'] = new_canonical_name

def recursive_rename(path, new_canonical_name, modify=False):
    sub_files = os.listdir(path)
    for sub_file in sub_files:
        target = os.path.join(path, sub_file)
        #rename directory and recurse down a level
        if os.path.isdir(target):
            old_name = sub_file.split('_')
            new_name = '_'.join([new_canonical_name] + old_name[1:])
            if modify:
                shutil.mv('{}'.format(os.path.join(path, sub_file), os.path.join(path, new_name)))
            else:
                print 'Moving:', os.path.join(path, sub_file), 'to', os.path.join(path, new_name)
            recursive_rename(target, new_canonical_name, modify=modify)
        #otherwise it is a file that needs to be renamed: two different procedures for
        #.xyz or for .json
        else:
            ext_str = os.path.splitext(sub_file)[1] 
            if ext_str == '.xyz':
                old_name = sub_file.split('_')
                new_name = '_'.join([new_canonical_name] + old_name[1:])
                if modify:
                    shutil.mv('{}'.format(os.path.join(path, sub_file), os.path.join(path, new_name)))
                else:
                    print 'Moving:', os.path.join(path, sub_file), 'to', os.path.join(path, new_name)
            elif ext_str == '.json':
                with open(os.path.join(path, 'subgb.json'),'r') as f:
                    subgb_dict = json.load(f)
                try:
                    subgb_dict['old_gbid'] = subgb_dict['gbid']
                except KeyError:
                    subgb_dict['old_gbid'] = subgb_dict['name']
                old_name = subgb_dict['old_gbid']
                subgb_dict['gbid'] = new_canonical_name +'_' + '_'.join(old_name[1:])
                if modify:
                    with open(os.path.join(path,'subgb.json')) as f:
                        json.dump(subgb_dict, f, indent=2)
                else:
                    print subgb_dict['gbid'], subgb_dict['old_gbid']
            elif ext_str == '.idx':
                pass
            elif ext_str == '.pbs':
                os.remove(os.path.join(path, sub_file)) 
            else:
                print 'unknown file type', path, sub_file, os.path.splitext(sub_file)

if __name__ == '__main__':
    #or_axes = ['001', '110', '111']
    #1
    #ENTER ORIENTATION AXIS LEVEL AND RENAME CANONICAL DIRS.
    #for or_x in or_axes:
    #    with pushd(or_x) as ctx0:
    #        rename_canonical_dirs()
    #2
    #Correct SubGrainBoundary Directories:
    or_axes = ['0_0_1', '1_1_0', '1_1_1']
    pot_dirs = ['POTBH','EAM_Dud','EAM_Men','EAM_Mish','EAM_Ack']
    for or_x in or_axes:
        with pushd(or_x) as ctx0:
            canonical_dirs = glob.glob('{}*'.format(or_x))
            canonical_dirs = filter(os.path.isdir, canonical_dirs)
            for canonical_dir in canonical_dirs:
                with pushd(canonical_dir) as ctx1:
                    print canonical_dir
                    with open('gb.json','r') as f:
                        gb_dict = json.load(f)
                    #use old canonical name:
                    #old_gbid =  gb_dict['old_gbid']
                    old_gbid =  gb_dict['gbid']
                    new_canonical_name = canonical_name_from_gb_json(gb_dict)
                    for pot in pot_dirs:
                        if os.path.isdir(pot):
                            with pushd(pot) as ctx2:
                                recursive_rename('./', new_canonical_name, modify=False)
                        else:
                           print 'No', pot, 'in', canonical_dir

