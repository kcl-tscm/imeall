import os
import glob
import pickle

class CrackDict(object):
    """
    Object to modify crack dictionaries. If from_name is true
    the property should be determined by some subset of characters from
    the directory name. Otherwise just specify a single name.
    The default is the temperature. Energy release rate
    would be another common one. Strain might also be a good idea.
    Params:
      `ps` starting char in job string
      `pf` finishing char in job string
      `pattern` job directory pattern
      `key`     key to add to dict
      `value`   value to associate with key
      `from_name` whether assign passed value or take value from job string
    """
    def __init__(self, key='sim_T', value=300.0, pattern="T*", from_name=True, ps=1, pf=4):
        self.pattern   = pattern
        self.from_name = from_name
        self.key       = key
        self.ps        = ps
        self.pf        = pf
        self.value     = 300.0

    def print_keys(self):
        """
        print keys of dicts in all directories matching
        pattern.
        """
        jobs = glob.glob(self.pattern)
        scratch = os.getcwd()
        for job in jobs:
            target = os.path.join(scratch, job)
            os.chdir(target)
            with open("crack_info.pckl",'r') as f:
                crack_dict = pickle.load(f)
            print job, crack_dict
        os.chdir(scratch)
        return

    def add_key_to_dict(self):
        """
        add key to dict in all directories matching pattern.
        """
        scratch = os.getcwd()
        jobs = glob.glob(self.pattern)
        for job in jobs:
            target = os.path.join(scratch, job)
            os.chdir(target)
            with open("crack_info.pckl",'r') as f:
                crack_dict = pickle.load(f)
            if self.from_name:
                print "adding key to ", job, "dict with value ", job[self.ps:self.pf]
                crack_dict[self.key] = float(job[self.ps:self.pf])
            elif not self.from_name:
                print "adding key to", job, "dict with value", self.value
                crack_dict[self.key] = self.value
            else:
                print "Something v. wrong"
                raise
            with open("crack_info.pckl",'w') as f:
                pickle.dump(crack_dict,f)
        os.chdir(scratch)
        return
