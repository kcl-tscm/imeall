import os
import sys
import json
import shutil
from   cStringIO       import StringIO
import ase.io        
from   ase.constraints import UnitCellFilter
from   ase.optimize    import BFGS, FIRE
from   quippy          import Atoms, Potential, frange, farray, fzeros
from   quippy.io       import AtomsWriter, AtomsReader, write
from   slabmaker.slabmaker  import build_tilt_sym_gb
import numpy as np
from   pprint import pprint
import argparse

class Capturing(list):
  def __enter__(self):
    self._stdout = sys.stdout
    sys.stdout=self._stringio = StringIO()
    return self

  def __exit__(self, *args):
    self.extend(self._stringio.getvalue().splitlines())
    sys.stdout = self._stdout

class ImeallIO(object):
  """
  Primary IO Object for routines for searching Imeall Directory tree
  creating new grains, and subgrain directories. 
  Sub-Grain directory (each 'massaged' grain boundary will have its own
  directory) supercells of the grain count as subgrains finally a subdirectory
  of defects is included for vacancies and interstitials.
  """
  def __init__(self):
#IO variables for VASP Calculations. 
#Might be better to have a separate VASP OBJECT templated POSCAR, INCAR Files
#then an Espresso Object Template
    self.vasp_template_dir =  '/projects/SiO2_Fracture/iron/vasp_template/'
    self.vasp_dict         = {'kpar'  :32, 'npar':16, 'magmom':3.0, 'n_at':'NOTSET','ediffg': -0.05}
    self.kpt_grid          = {'kx':12, 'ky':12, 'kz':1}
    self.runsh             = {'nodes':512, 'time':360}

  def make_dir(self, target_dir, dir_name):
    target_subdir = os.path.join(target_dir, dir_name)
    if not os.path.isdir(target_subdir):
      os.makedirs(target_subdir)
      print 'Created {0}'.format(target_subdir)
    else:
      print '\t directory already exists'
    return target_subdir

  def load_json(self, json_file):
		with open(json_file, 'r') as datfile:
			gb_data = json.load(datfile)
		return gb_data

  def find_subdir(self, target_dir, calc_suffix):
    for dir in os.listdir(target_dir):
      if calc_suffix in dir:
        from_dir = os.path.join(target_dir, dir)
        from_dir_name = dir
        return from_dir, from_dir_name
    print 'no directory with calc_suffix,', calc_suffix
    return None
# Common pattern might be I have a list of grain dir directories
# with xyz files in them. Say the relaxed EAM grainboundaries for
# a certain orientation. /gb/EAM/gbid_criteria/gbid_traj.xyz.


  def copy_struct(self, from_dir, target_dir, from_sub_dir, 
      target_sub_dir, calc_suffix='_d2.0', calc_point='final', md_step=1):
# The calc_suffix determines which subgrain directories to
# copy across. For instance _d2.0 will search for all subgrains
# generated with the deletion critera for atoms with nn distance < 2.0 A.
# This is useful if there are a number of subgrains with different deletion
# criteria or defect concentrations etc. calculated in the EAM directory
# and we want to copy those over to a DFT calc, or a Tightbinding representation of the
# grain boundary. Another intended use is if I have a
# eam grain boundary or cleavage plane
# md run and I want to copy across snap 
# shots of the trajectory to do some DFT on.
    json_file = os.path.join(from_dir, 'gb.json')
    gb_data = self.load_json(json_file)
    from_dir = os.path.join(from_dir, from_sub_dir)	
    calc_dirs = os.listdir(from_dir)
#Search for path that matches the calc_suffix
    for dir in calc_dirs:
      if calc_suffix in dir:
        from_dir = os.path.join(from_dir, dir)
        from_dir_name = dir
        break
    if not os.path.isdir(from_dir):
      print 'no directory with calc_suffix,', calc_suffix
# Can start grain_boundary at arbitrary point in the trajectory
# start at the end of the md run, midpoint, beginning, or arbitrary, 
# position 
    at = AtomsReader('{0}/{1}_traj.xyz'.format(from_dir, gb_data['gbid']))
    if calc_point =='final':
      grain = at[-1]
    elif calc_point =='midpoint':
      grain = at[int(float(len(at)/2.))]
    elif calc_point =='beginning':	
      grain = at[1]
    elif calc_point =='arbitrary':	
      grain = at[md_step]
    else:
      print 'invalid starting point'
# Now write atom object to target directory:
    target_dir = os.path.join(target_dir, os.path.join(target_sub_dir, from_dir_name))
    if not os.path.isdir(target_dir):
        os.makedirs(target_dir)
    else:
      print target_dir, ' already exists.'
# We have now copied the desired xyz (last snap shot of trajector)
# file from from_dir/EAM
# to target_dir/DFT along with the gb json data file.
    out = AtomsWriter('{0}'.format(os.path.join(target_dir,
                      '{0}.xyz'.format(gb_data['gbid']))))
    out.write(grain)
    gb_data['n_at'] = len(grain)
    with open(os.path.join(target_dir,'gb.json'),'w') as json_file:
      json.dump(gb_data, json_file, indent=2)

  def xyz_to_vasp(self, target_dir, target_sub_dir, calc_suffix='_d2.0'):
# Initialize a vasp calculation. Read in the xyz file 
# if it is present and generate an INCAR, POSCAR, 
# run_vasp file in the Directory the vasp template 
    json_file = os.path.join(target_dir, 'gb.json')
    gb_data = self.load_json(json_file)
    target_dir = os.path.join(target_dir, target_sub_dir)
    calc_dirs = os.listdir(target_dir)
    for dir in calc_dirs:
      if calc_suffix in dir:
        target_dir = os.path.join(target_dir, dir)
        break
    if not os.path.isdir(target_dir):
      print 'no directory with calc_suffix,', calc_suffix
    incar_tmp_file = os.path.join(self.vasp_template_dir, 'INCAR')
    with open(incar_tmp_file, 'r') as incar_file:
      incar_template = incar_file.read()

    grain = Atoms('{0}'.format(os.path.join(target_dir,
    '{0}.xyz'.format(gb_data['gbid']))))
    n_at           = len(grain)
#Write INCAR file to target directory
    self.vasp_dict['n_at'] = n_at
    with open(os.path.join(target_dir, 'INCAR'),'w') as incar_file:
      print >> incar_file, incar_template.format(**self.vasp_dict)
#Write run.sh with job time and number of nodes to job directory
    runsh_tmp_file = os.path.join(self.vasp_template_dir, 'run.sh')
    with open(runsh_tmp_file, 'r') as runsh_file:
      runsh_template = runsh_file.read()
    with open(os.path.join(target_dir, 'run.sh'),'w') as runsh_file:
      print >> runsh_file, runsh_template.format(**self.runsh)
      os.chmod(os.path.join(target_dir, 'run.sh'), 0775)
#Write  POSCAR to target directory
    ase.io.write(os.path.join(target_dir, 'POSCAR'), grain)	
#Write KPOINT file to target directory
    with open(os.path.join(self.vasp_template_dir, 'KPOINTS'), 'r') as kpt_tmp_file:
      kpt_template = kpt_tmp_file.read() 
    with open(os.path.join(target_dir, 'KPOINTS'), 'w') as kpt_file:
      print >>kpt_file, kpt_template.format(**self.kpt_grid)
#Copy run_vasp configuration file and the pseudopotential
    shutil.copy(os.path.join(self.vasp_template_dir, 'run_vasp'), target_dir)
    shutil.copy(os.path.join(self.vasp_template_dir, 'POTCAR'), target_dir)

  def submit_job(self, target_dir, target_sub_dir, calc_suffix='_d2.0'):
    target_dir = os.path.join(target_dir, target_sub_dir)
#Identify the subdirectory by the calculation suffix:
    target_dir, target_dir_name = self.find_subdir(target_dir, calc_suffix)
    try:
      gb_data = json.load(open(os.path.join(target_dir,'gb.json')))
      pprint(gb_data)
    except:
      print 'Error Opening JSON FILE'
    try:
#os.system('cd {0}; ./run.sh'.format(target_dir))
      pass
    except:
      print 'Error submitting ', target_dir

class GBRelax(object):
  def __init__(self, grain_dir='./', gbid='0000000000', calc_type='EAM',
               potential = 'IP EAM_ErcolAd', param_file = 'iron_mish.xml',
               traj_file='traj.xyz'):
# Here we initialize some naming conventions, the calculation type, and
# necessary input files.
    self.gbid        =  gbid
# grain_dir is the overarching grain directory (the theme):
# the target dir is the subdirectory depending on the
# calculation type (flavour) at the highest level and then the
# variations (deletions, translation, substitutions, vacancies) 
# on the theme are subgrain(s) i.e. in subgrain_dir named
# according to the variation. atom deletion for atoms within radius
# rcut is written gbid_r2.0, translations gbid_tx_0.1, gbid_ty,0.2 etc.
# Hydrogen occlusion is denoted interms of concentration, or number of H?
# That remains to be seend
    self.grain_dir    = grain_dir
    self.calc_dir   = os.path.join(grain_dir, calc_type)
    self.subgrain_dir =''
    self.name        = 'subgrain'
    self.calc_type   = calc_type
    self.struct_file = ''
    self.param_file  = param_file
    self.potential   = potential
    self.fmax        = 0.5E-2
    self.traj_file   = traj_file

    print 'Generate Job in:'
    print self.grain_dir
    print self.calc_dir

  def gen_super_rbt(self, bp=[],v=[], rbt=[0.0, 0.0], sup_v=6, sup_bxv=4, rcut=2.0):
    '''
      Create a grain boundary with rigid body translations.
    '''
    io = ImeallIO()
    grain = build_tilt_sym_gb(bp=bp, v=v, rbt=rbt)
# For RBT we build a top level dir with just the translated supercell and no
# deleted atoms then we create subdirectories with particular deletion criteria
# to them. This means we don't have to duplicate the structure files.
    struct_dir = os.path.join(self.grain_dir,'structs')
    struct_files = os.listdir(struct_dir)
    loc_name    = '{0}_v{1}bxv{2}_tv{3}bxv{4}_d{5}z'.format(self.gbid,
    str(sup_v), str(sup_bxv), str(round(rbt[0],2)), str(round(rbt[1],2)), str(rcut))
#Check if this initial grain with the atom deletion criterion and the rigid body translations
#is present if not generate it:
    print loc_name+'.xyz'
    #if loc_name+'.xyz' not in struct_files:
    #Quite a few non-supercells were in the directories so regenerating from scratch.
    if True:
      print 'GENERATING STRUCTURE'
      tmp0, tmp1, grain  = self.gen_super(grain=grain, rbt=rbt, sup_v=sup_v, sup_bxv=sup_bxv, rcut=0.0)
      self.name    = '{0}_v{1}bxv{2}_tv{3}bxv{4}'.format(self.gbid,
      str(sup_v),str(sup_bxv), str(round(rbt[0],2)), str(round(rbt[1],2)))
      self.subgrain_dir = io.make_dir(self.calc_dir, self.name)
      grain.write('{0}.xyz'.format(os.path.join(self.subgrain_dir, self.name)))
# Now create sub_dir with the appropriate deletion criteria:
      self.name    = '{0}_v{1}bxv{2}_tv{3}bxv{4}_d{5}z'.format(self.gbid,
      str(sup_v), str(sup_bxv), str(round(rbt[0],2)), str(round(rbt[1],2)), str(rcut))
      self.subgrain_dir = io.make_dir(self.subgrain_dir, self.name)
      grain = self.delete_atoms(grain=grain, rcut=rcut)
#Deposit all initial structures in the struct dir:
      grain.write('{0}.xyz'.format(os.path.join(struct_dir, self.name)))
    else:
      print 'STRUCTURE ALREADY GENERATED'
      self.name    = '{0}_v{1}bxv{2}_tv{3}bxv{4}'.format(self.gbid,
      str(sup_v),str(sup_bxv), str(round(rbt[0],2)), str(round(rbt[1],2)))
      self.subgrain_dir = io.make_dir(self.calc_dir, self.name)
      self.name    = '{0}_v{1}bxv{2}_tv{3}bxv{4}_d{5}z'.format(self.gbid,
      str(sup_v), str(sup_bxv), str(round(rbt[0],2)), str(round(rbt[1],2)), str(rcut))
      self.subgrain_dir = io.make_dir(self.subgrain_dir, self.name)

# Finally deposit json and grain file with translation information.
    try:
      f = open('{0}/subgb.json'.format(self.subgrain_dir), 'r')
      j_dict = json.load(f)
      f.close()
    except IOError:
       f = open('{0}/subgb.json'.format(self.subgrain_dir), 'w')
       j_dict = {}
#terms to append to subgrain dictionary:
    j_dict['param_file'] = self.param_file
    j_dict['name'] = self.name
    j_dict['rbt']  = rbt
    j_dict['rcut'] = rcut

    f = open('{0}/subgb.json'.format(self.subgrain_dir), 'w')
    json.dump(j_dict,f, indent=2)
    f.close()


  def gen_super(self, grain=None, rbt=None, sup_v=6, sup_bxv=2, rcut=2.0):
    ''' 
       To create a grain boundary super cell we use the parameters of
       Rittner and Seidman (PRB 54 6999).
    '''
    io = ImeallIO()
    if rbt == None:
      x = Atoms('{0}.xyz'.format(os.path.join(self.grain_dir, self.gbid)))
    else:
      x = Atoms(grain)

    struct_dir = os.path.join(self.grain_dir, 'structs')  
    self.name  = '{0}_v{1}bxv{2}_tv{3}bxv{4}_d{5}z'.format(self.gbid,
    str(sup_v), str(sup_bxv), '0.0', '0.0', str(rcut))
    struct_files = os.listdir(struct_dir)
# If this structure already exists for this grain boundary
# don't recreate it. Otherwise build the supercell and deposit 
# it in the structs directory.
    #if self.name+'.xyz' in struct_files:
    if False:
      print 'Structure file already exists'
      pass
    else:
      if rcut > 0.0:
        x.set_cutoff(2.4)
        x.calc_connect()
        x.calc_dists()
        rem=[]
        r = farray(0.0)
        u = fzeros(3)
        for i in frange(x.n):
          for n in frange(x.n_neighbours(i)):
            j = x.neighbour(i, n, distance=3.0, diff=u)
            if x.distance_min_image(i,j) < rcut and j!=i:
              rem.append(sorted([j,i]))
        rem = list(set([a[0] for a in rem]))
        if len(rem) >0:
          x.remove_atoms(rem)
        else:
          print 'No duplicate atoms in list.'
      else:
        pass
  # Now create super cell:
      x = x*(sup_v, sup_bxv, 1)
      x.set_scaled_positions(x.get_scaled_positions())
      x.write('{0}.xyz'.format(os.path.join(struct_dir, self.name)))

    if rbt == None:
      self.struct_file  = self.name
      self.subgrain_dir = io.make_dir(self.calc_dir, self.name)
      try:
        with  open('{0}/subgb.json'.format(self.subgrain_dir), 'r') as f:
          j_dict = json.load(f)
      except IOError:
        j_dict             = {}
      j_dict['name']       = self.name
      j_dict['param_file'] = self.param_file
      j_dict['rbt']        = [0.0, 0.0]
      j_dict['rcut']       = rcut
      with open('{0}/subgb.json'.format(self.subgrain_dir), 'w') as f:
        json.dump(j_dict, f, indent=2)
#Deposit structure in the structs directory of the grain
    else:
      return sup_v, sup_bxv, x

  def delete_atoms(self, grain=None, rcut=2.0):
# Delete atoms below a certain threshold
    io = ImeallIO()
    if grain == None:
      x = Atoms('{0}.xyz'.format(os.path.join(self.grain_dir, self.gbid)))
    else:
      x = Atoms(grain)
    x.set_cutoff(2.4)
    x.calc_connect()
    x.calc_dists()
    rem=[]
    r = farray(0.0)
    u = fzeros(3)
    for i in frange(x.n):
      for n in frange(x.n_neighbours(i)):
        j = x.neighbour(i, n, distance=3.0, diff=u)
        if x.distance_min_image(i, j) < rcut and j!=i:
          rem.append(sorted([j,i]))
    rem = list(set([a[0] for a in rem]))
    if len(rem) > 0:
      x.remove_atoms(rem)
    else:
      print 'No duplicate atoms in list.'
    if grain == None:
      self.name    = '{0}_d{1}'.format(self.gbid, str(rcut)) 
      self.subgrain_dir = io.make_dir(self.calc_dir, self.name)
      self.struct_file  = gbid + '_' + 'n' + str(len(rem)) + 'd' + str(rcut)  
      x.write('{0}.xyz'.format(os.path.join(self.subgrain_dir, self.struct_file)))
      return len(rem)
    else:
      return x

  def gen_pbs(self, time='02:30:00', queue='serial.q'):
    ''' 
      Generates a pbs file from the template stored in pbs_str. Writes 
    '''
    pbs_str = open('/users/k1511981/pymodules/templates/calc_ada.pbs','r').read()
    pbs_str = pbs_str.format(jname='fe'+self.name, xyz_file='{0}.xyz'.format(self.name), 
                             time=time, queue=queue)
    print os.path.join(self.subgrain_dir, 'fe{0}.pbs'.format(self.name))
    with open(os.path.join(self.subgrain_dir, 'fe{0}.pbs'.format(self.name)) ,'w') as pbs_file:
      print >> pbs_file, pbs_str

  def translate(self):
    pass
# Displace grain_a relative to grain_b
# can be accomplished by moving along x or y
# and wrapping atoms back to unit cell
# Take the grain_boundary structure and relax it
  def go_relax(self, maxiters = 400):
    if self.struct_file == '':
      print 'No pointer to struct file. Have you specified a job?'
      return None
    else:
      print '\t Print Relaxing {0}.xyz with {1} potential'.format(self.struct_file, self.calc_type)
      print '\t In directory {0}'.format(self.subgrain_dir) 
    if self.calc_type[:3] == 'EAM':
      grain       = Atoms('{0}.xyz'.format(os.path.join(self.subgrain_dir,
                            self.struct_file)))
# Load potential and Attach calculator
      pot         = Potential(self.potential, param_filename=self.param_file)
      grain.set_calculator(pot)
      strain_mask = [0,0,0,0,0,0]
      ucf            = UnitCellFilter(grain, mask=strain_mask)
      self.traj_file = self.name+'.traj'
      traj_loc    = os.path.join(self.subgrain_dir, self.traj_file)
      traj        = ase.io.Trajectory(traj_loc, 'w', grain)
      E_gb_init   = grain.get_potential_energy()
      opt         = BFGS(ucf)
      opt.attach(traj.write, interval=25)
#Cant get printing the output from bfgs to behave...
      #with open(os.path.join(self.subgrain_dir, 'bfgsmin.txt'), 'w') as outfile:
      #  sysold = sys.stdout
      #  sys.stdout = StringIO()
      opt.run(fmax = self.fmax, steps=maxiters)
      #  print >> outfile, output
      #  sys.stdout.close()
      #  sys.stdout = sysold
      traj.close()
      out = AtomsWriter('{0}'.format(os.path.join(self.subgrain_dir,
                        '{0}_traj.xyz'.format(self.gbid))))
#This is a hack way of making sure we only have xyz files of trajectories
      for at in AtomsReader(traj_loc):
        out.write(at)
      out.close()
      os.remove(traj_loc)
#convert traj file to extended xyz
      E_gb = grain.get_potential_energy()
      cell = grain.get_cell()
      A    = cell[0][0]*cell[1][1]
# Calculation dumps total energyenergy and grainboundary area data to json file.
      gb_dict = {'E_gb':E_gb, 'E_gb_init':E_gb_init, 'A': A, 'n_at':len(grain)}
      with open(os.path.join(self.subgrain_dir, 'subgb.json'), 'w') as outfile:
        json.dump(gb_dict, outfile, indent=2)

if __name__=='__main__':
# string for orientation axis e.g. '110'
# run_dyn parsing builds a command line tool for managing grain boundary
# database
  parser = argparse.ArgumentParser()
  parser.add_argument("-p", "--prefix", help = "Subsequent commands will act on all \
                                                subdirectories with first characters matching prefix.", required=True)
  parser.add_argument("-ct", "--calc_type", help = "Name of calculation type TB, EAM, DFT, etc.", required=True)
  parser.add_argument("-q",  "--queue", help = "Jobs will be submitted to this queue.", default='smp.q')
  parser.add_argument("-t",  "--time",  help = "Time limit on jobs.", default='1:00:00')
  parser.add_argument("-rc", "--rcut", type=float, help = "Deletion criterion for nearest neighbour atoms.", default=2.0)

  args = parser.parse_args()
  calc_type = args.calc_type
  prefix    = args.prefix
  queue     = args.queue
  time      = args.time
  rcut      = args.rcut
  print rcut, type(rcut)

# Each calculation type is associated with a potential:
  if calc_type == 'EAM_Mish':
    param_file = 'iron_mish.xml'
  elif calc_type == 'EAM_Men':
    param_file = 'Fe_Mendelev.xml'
  elif calc_type == 'EAM_Ack':
    param_file = 'Fe_Ackland.xml'
  elif calc_type == 'EAM_Dud':
    param_file = 'Fe_Dudarev.xml'
  else:
    print 'No available potential corresponds to this calculation type.'
    sys.exit()

  jobdirs = []
  for thing in os.listdir('./'):
    if os.path.isdir(thing) and thing[:len(prefix)]==prefix:
      jobdirs.append(thing)

  for job_dir in jobdirs[:]:
    gbid    = job_dir.strip('/')
    print '\n'
    print '\t', gbid
    print '\n'
    gbrelax = GBRelax(grain_dir=job_dir, gbid=gbid, calc_type=calc_type,
                      potential = 'IP EAM_ErcolAd', param_file=param_file)
    sup_v   = 6
    sup_bxv = 2
#check if structure already exists in structs file:
    with open(os.path.join(job_dir,'gb.json')) as f:
      grain_dict = json.load(f)
    bp = grain_dict['boundary_plane']
    v  = grain_dict['orientation_axis']
    for i in np.linspace(0.25,1.0, 4):
      for j in np.linspace(0.125,1.0,8):
        gbrelax.gen_super_rbt(rcut=rcut, bp=bp, v=v, sup_v=sup_v, sup_bxv=sup_bxv, rbt=[i,j])
        gbrelax.gen_pbs(time=time, queue=queue)

#######################################
#######################################
## RELAX EAM FOR LIST OF DIRECTORIES ##
#######################################
#######################################
#  for job_dir in jobdirs[:]:
#    gbid    = job_dir.strip('/')
#    print '\n'
#    print '\t', gbid
#    print '\n'
#    gbrelax = GBRelax(grain_dir=job_dir, gbid=gbid, calc_type='EAM', 
#                      potential = 'IP EAM_ErcolAd', param_file='./iron_mish.xml')
#    gbrelax.delete_atoms()
#    gbrelax.go_relax()
#########################################################################
#########################################################################
#print jobdirs
#for job_dir in jobdirs[:]:
#  print job_dir
#  gbid    = job_dir.strip('/')
#  gbrelax = GBRelax(grain_dir=job_dir, gbid=gbid, calc_type='EAM_Men', 
#                    potential='IP EAM_ErcolAd', param_file='Fe_Mendelev.xml')
#                    potential = 'IP EAM_ErcolAd', param_file='./iron_mish.xml')
#  with open(os.path.join(job_dir,'gb.json')) as f:
#    grain_dict = json.load(f)
#  bp = grain_dict['boundary_plane']
#  v  = grain_dict['orientation_axis']
#  for i in np.linspace(0.25,1.0,4):
#    for j in np.linspace(0.125,1.0,8):
#      gbrelax.gen_super_rbt(bp=bp, v=v, rbt=[i,j], rcut=2.0)
#      gbrelax.gen_pbs()
#########################################################################
#########################################################################
## COPY Directories across im_io.copy_struct(dir, sub_dir, dir, sub_dir)#
#########################################################################
#########################################################################
##  im_io = ImeallIO()
##  for job_dir in jobdirs:
##    print job_dir
##    im_io.copy_struct(job_dir, job_dir, 'EAM', 'DFT')
#########################################################################
#########################################################################
## GENERATE list of poscar files etc in the sub_directories            ##
#########################################################################
#########################################################################
#  im_io = ImeallIO()
#  for job_dir in jobdirs[1:]:
#    print job_dir
#    im_io.xyz_to_vasp(job_dir,'DFT')
############################################################################
############################################################################
## Our final pattern is to scan through all the subgrains and resubmit    ##
## unless a particular convergence has been achieved.                     ##
## This is going to require a few intelligent parsers of the VASP OSZICAR ##
## And an autoroutine for                                                 ##
############################################################################
#  im_io = ImeallIO()
#  for job_dir in jobdirs[1:]:
#    print job_dir
#    im_io.submit_job(job_dir, 'DFT')

