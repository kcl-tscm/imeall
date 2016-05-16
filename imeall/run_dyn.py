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

from pprint import pprint

class Capturing(list):
  def __enter__(self):
    self._stdout = sys.stdout
    sys.stdout=self._stringio = StringIO()
    return self

  def __exit__(self, *args):
    self.extend(self._stringio.getvalue().splitlines())
    sys.stdout = self._stdout

class ImeallIO(object):
# Primary IO Object for routines for searching Imeall Directory tree
# creating new grains, and subgrain directories. 
# Sub-Grain directory (each 'massaged' grain boundary will have its own
# directory) supercells of the grain count as subgrains finally a subdirectory
# of defects is included for vacancies and interstitials.
  def __init__(self):
#IO variables for VASP Calculations. 
#Might be better to have a separate VASP OBJECT templated POSCAR, INCAR Files
#then an Espresso Object Template
    self.vasp_template_dir = '/projects/SiO2_Fracture/iron/vasp_template/'
    self.vasp_dict         = {'kpar'  :32, 'npar':16, 'magmom':3.0, 'n_at':'NOTSET',
					      		        	'ediffg': -0.05}
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
# copy across. For instance _d2.0 which search for all subgrains
# generated with the deletion critera for atoms of anything under 2.0 A.
# This is useful if there are a number of subgrains with different deletion
# criteria or defect concentrations etc. calculated in the EAM directory
# and we want to copy those over to a DFT calc, or a Tightbinding representation of the
# grain boundary. Another intended use is if I have a eam grain boundary or cleavage plane
# md run and I want to copy across snap shots of the trajectory to do some DFT on.
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

  def gen_super(self, rcut=2.0):
# Create a grainboundary super cell we use the parameters of
# Rittner and Seidman.
    io = ImeallIO()
    x = Atoms('{0}.xyz'.format(os.path.join(self.grain_dir, self.gbid)))
    x.set_cutoff(3.0)
    x.calc_connect()
    x.calc_dists()
    rem=[]
    r = farray(0.0)
    u = fzeros(3)
    for i in frange(x.n):
      for n in frange(x.n_neighbours(i)):
      	j = x.neighbour(i, n, distance=3.0, diff=u)
        if 0. < x.distance_min_image(i,j)< rcut and j!=i:
       	  rem.append(sorted([j,i]))
    rem = list(set([a[0] for a in rem]))
    if len(rem) >0:
      x.remove_atoms(rem)
    else:
      print 'No duplicate atoms in list.'
# Now create super cell
    n = 2
    m = 6
    x = x*(m,n,1)
# Now Generate Subgrain directory
    x.set_scaled_positions(x.get_scaled_positions())
    self.name    = '{0}_v{1}bxv{2}_d{3}z'.format(self.gbid, str(m),str(n),str(rcut)) 
    self.struct_file  = self.name
    self.subgrain_dir = io.make_dir(self.calc_dir, self.name)
    x.write('{0}.xyz'.format(os.path.join(self.subgrain_dir, self.name)))

  def delete_atoms(self, rcut=2.0):
# Delete atoms below a certain threshold
    io = ImeallIO()
    x = Atoms('{0}.xyz'.format(os.path.join(self.grain_dir, self.gbid)))
    x.set_cutoff(3.0)
    x.calc_connect()
    x.calc_dists()
    rem=[]
    for j in frange(x.n):
      for i in frange(j, x.n):
        if 0. < x.distance_min_image(i, j) < rcut and j!=i:
          rem.append(sorted([j,i]))
    rem = list(set([a[0] for a in rem]))
    if len(rem) >0:
      x.remove_atoms(rem)
    else:
      print 'No duplicate atoms in list.'
    self.name    = '{0}_d{1}'.format(self.gbid, str(rcut)) 
    self.subgrain_dir = io.make_dir(self.calc_dir, self.name)
    self.struct_file  = gbid + '_' + 'n' + str(len(rem)) + 'd' + str(rcut)  
    x.write('{0}.xyz'.format(os.path.join(self.subgrain_dir, self.struct_file)))
    return len(rem)

  def gen_pbs(self, time='12:00:00'):
		pbs_str = open('/users/k1511981/pymodules/templates/calc_ada.pbs','r').read()
		pbs_str = pbs_str.format(jname = 'fe'+self.name, xyz_file='{0}.xyz'.format(self.struct_file), time=time)
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
    if self.calc_type == 'EAM':
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
  jobdirs = []
#string for orientation axis e.g. '110'
  or_string = sys.argv[1]
  for thing in os.listdir('./'):
    if os.path.isdir(thing) and thing[:3]==or_string:
      print thing
      jobdirs.append(thing)
#######################################
#######################################
## RELAX EAM FOR LIST OF DIRECTORIES ##
#######################################
#######################################
##  for job_dir in jobdirs[:]:
##    gbid    = job_dir.strip('/')
##    print '\n'
##    print '\t', gbid
##    print '\n'
##    gbrelax = GBRelax(grain_dir=job_dir, gbid=gbid, calc_type='EAM', 
##                      potential = 'IP EAM_ErcolAd', param_file='./iron_mish.xml')
##    gbrelax.delete_atoms()
##    gbrelax.go_relax()
#########################################################################
#########################################################################
print jobdirs
for job_dir in jobdirs[:]:
  gbid    = job_dir.strip('/')
  print '\n'
  print '\t', gbid
  print '\n'
  gbrelax = GBRelax(grain_dir=job_dir, gbid=gbid, calc_type='EAM', 
                    potential = 'IP EAM_ErcolAd', param_file='./iron_mish.xml')
  gbrelax.gen_super(rcut=2.3)
  gbrelax.gen_pbs()
#  gbrelax.go_relax()
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

