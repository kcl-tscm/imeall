import argparse
import subprocess
import datetime

class SyncDB(object):
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
#Use rsync to sync the db with that server
    #sync_arg = "rsync -auv --exclude '*/Fracture/*' --exclude-from 'rsync_exclude.txt' lambert@mira.alcf.anl.gov:/home/lambert/iron/grain_boundaries ./"
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
  parser.add_argument("-s", "--server", help="name of server to synce with")

  mira_params = dict(sync_log="db_synclog", exclude="'*/Fracture/*'", exclude_from="rsync_exclude.txt", rsync_args="-auv",
                     src="lambert@mira.alcf.anl.gov:/home/lambert/iron/grain_boundaries", target="./") 

  sync = SyncDB(**mira_params)
  sync.sync_db(server="mira")




