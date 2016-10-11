import datetime
import argparse
import subprocess

class SyncDB(object):
  """
    Class for synchronizing the local database with remote servers e.g. ADA, Rosalind, Archer, MIRA etc.
    Needs to be updated so that subprocess can handle log in credentials via the browser.
    Database information is stored in db_synclog.
  """
  def __init__(self, sync_log="db_synclog",exclude="'*/Fracture/*'", exclude_from="rsync_exclude.txt",
               src="lambert@mira.alcf.anl.gov:/home/lambert/iron/grain_boundaries",
               rsync_args="-auv", target="./"):

    self.src          = src
    self.target       = target
    self.exclude      = exclude
    self.exclude_from = exclude_from
    self.sync_log     = sync_log
    self.rsync_args   = rsync_args

  def sync_db(self, server="mira"):
# Use rsync to sync the db with that server
# sync_arg = "rsync -auv --exclude '*/Fracture/*' --exclude-from 'rsync_exclude.txt' 
# lambert@mira.alcf.anl.gov:/home/lambert/iron/grain_boundaries ./"
    sync_args = ["rsync"] + [self.rsync_args] + ["--exclude"] + [self.exclude] + ["--exclude-from"] + [self.exclude_from] + [self.src] + [self.target]
    print sync_args
    with open(self.sync_log,"a") as output:
      right_now = datetime.datetime.today().strftime("%H:%M:%S %d-%m-%Y")
      output.write("\n\n")
      output.write("Syncing local db with {server} ".format(server=server) + right_now +"\n" )
      output.flush()
      job = subprocess.Popen(sync_args, stdout=output)
      job.wait()

if __name__=='__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument("-s", "--server", help="name of server to sync with")
  parser.add_argument("-a", "--ada",  action="store_true")
  parser.add_argument("-m", "--mira", action="store_true")
  args = parser.parse_args()

  mira_params = dict(sync_log="db_synclog", exclude="'*/Fracture/*'", exclude_from="rsync_exclude.txt", rsync_args="-av",
                     src="lambert@mira.alcf.anl.gov:/home/lambert/iron/grain_boundaries", target="./") 

  ada_params = dict(sync_log="db_synclog", exclude="'*/Fracture/*'", exclude_from="rsync_exclude.txt", rsync_args="-auv", 
                     src="k1511981@ada.hpc.kcl.ac.uk:/users/k1511981/sharedscratch/grain_boundaries/grain_boundaries", target="./")

  #To sync the directory structure of the grainboundary
  #rsync -a -f"+ */" -f"- *" source/ destination/
  #rsync_args='-a -f"+ */" -f"- *'.split()
  #http://superuser.com/questions/156664/what-are-the-differences-between-the-rsync-delete-options
  emptydir_params = dict(sync_log="db_synclog", exclude="'*/Fracture/*'", exclude_from="rsync_exclude.txt", rsync_args="-auv", src="", target="")

  if args.mira:
    sync_mira = SyncDB(**mira_params)
    sync_mira.sync_db(server="mira")
  
  if args.ada:
    sync_ada  = SyncDB(**ada_params)
    sync_ada.sync_db(server="ada")




