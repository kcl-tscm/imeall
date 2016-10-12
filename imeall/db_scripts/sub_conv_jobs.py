import os
import json
import subprocess
from imeall.models import GBAnalysis

analysis = GBAnalysis()
converged, unconverged, missing_key = analysis.list_unconverged()

print 'Unconverged'
scratch = os.getcwd()
for gb in unconverged[5:]:
  print os.path.join(gb[0], gb[1])
# Typical pattern for submitting jobs:
# changedir_to_jobdir -> determine_params_from_json -> modify_pbs_file -> submit job
  job_dir = os.path.join(gb[0], gb[1])
  os.chdir(os.path.join(scratch, job_dir))
  with open('subgb.json','r') as f:
    gb_dict = json.load(f)
# pbs_str   = open('/users/k1511981/pymodules/templates/calc_rundyn.pbs', 'r').read()
# gb_args   = '-rc {rc} -i_v {i_v} -i_bxv {i_bxv}'.format(rc=gb_dict['rcut'], i_v=gb_dict['rbt'][0], i_bxv=gb_dict['rbt'][1])
  pbs_str    = open('/users/k1511981/pymodules/templates/relax.pbs', 'r').read()
  input_file = gb_dict['gbid']
  pbs_str = pbs_str.format(jname='fe'+os.path.basename(gb[1][:8]), time='1:00:00', queue='smp.q', input_file=input_file)
  with open('relax.pbs', 'w') as pbs_file:
    print >> pbs_file, pbs_str
  qsub_args = 'qsub relax.pbs'
  job = subprocess.Popen(qsub_args.split())
  job.wait()
