import datetime
import argparse
import subprocess

class SyncDB(object):
  """
  :class:`SyncDB` for synchronizing the local database with remote servers.
  Updated file information is stored in db_synclog text file.
  Initialize object with all parameters required to perform the synchronization.
  Local modifications of this file to add relevant servers.
  e.g.
    sync_arg = "rsync -auv --exclude '*/Fracture/*' --exclude-from 'rsync_exclude.txt' 
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
  """
  :method:`sync_db` Use rsync to sync the db with that server.
  """
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
  parser.add_argument("-n", "--dryrun", help="perform a dry run.", action="store_true")
  parser.add_argument("-s", "--server", help="name of server to sync with.")
  parser.add_argument("-a", "--ada",    help="sync with ada.",  action="store_true")
  parser.add_argument("-m", "--mira",   help="sync with mira.", action="store_true")
  parser.add_argument("-r", "--rosa",   help="sync with rosalind.", action="store_true")
  #parser.add_argument("-ns", "--new_server", help="sync with new_server.", action="store_true")
  parser.add_argument("-f", "--fracrosa",   help="sync only fracture files from rosalind.", action="store_true")
  parser.add_argument("-e", "--exclude", help="exclude the following file.")
  args = parser.parse_args()

  rsync_args = "-auv"
  if args.dryrun:
    rsync_args +='n'

#Example synchronization dictionaries. Follow this pattern to create unique synchronzation dictionaries:
  #newserver_params = dict(sync_log="db_synclog", exclude="'*/Fracture/*'", exclude_from="rsync_exclude.txt", rsync_args="-av",
  #                   src="new_user@personal_server:/grain_boundaries_db", target="./") 
  #if args.new_server:
  #  sync_mira = SyncDB(**mira_params)
  #  sync_mira.sync_db(server="mira")
  mira_params = dict(sync_log="db_synclog", exclude="'*/Fracture/*'", exclude_from="rsync_exclude.txt", rsync_args="-av",
                     src="lambert@mira.alcf.anl.gov:/home/lambert/iron/grain_boundaries", target="./") 

  ada_params  = dict(sync_log="db_synclog", exclude="'*/Fracture/*'", exclude_from="rsync_exclude.txt", rsync_args="-auv", 
                     src="k1511981@ada.hpc.kcl.ac.uk:/users/k1511981/sharedscratch/grain_boundaries/grain_boundaries", target="./")

  rosa_params = dict(sync_log="db_synclog", exclude="'*/Fracture/*'", exclude_from="rsync_exclude.txt", rsync_args="-auv", 
                     src="k1511981@login.rosalind.compute.estate:/users/k1511981/sharedscratch/grain_boundaries", target="./")

  frac_params = dict(sync_log="db_synclog", exclude="'*/Fracture/*'", exclude_from="rsync_exclude.txt", rsync_args="-auv",
                     src="k1511981@login.rosalind.compute.estate:/users/k1511981/sharedscratch/grain_boundaries/alphaFe/f*",
                     target="./grain_boundaries/alphaFe/000/")


  emptydir_params = dict(sync_log="db_synclog", exclude="'*/Fracture/*'", exclude_from="rsync_exclude.txt", rsync_args="-auv", src="", target="")

  if args.mira:
    sync_mira = SyncDB(**mira_params)
    sync_mira.sync_db(server="mira")
  
  if args.ada:
    sync_ada  = SyncDB(**ada_params)
    sync_ada.sync_db(server="ada")

  if args.rosa:
    sync_ada  = SyncDB(**rosa_params)
    sync_ada.sync_db(server="rosalind")

  if args.fracrosa:
    sync_ada  = SyncDB(**frac_params)
    sync_ada.sync_db(server="rosalind")
