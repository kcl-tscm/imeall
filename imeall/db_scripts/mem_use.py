import os
import numpy as np

def get_size(start_path = '.'):
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(start_path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            total_size += os.path.getsize(fp)
    return total_size


print '\t Total memory usage in this and nested directories:'
print '\t', round(float(get_size())/float(np.power(2,30)),2), 'Gb'
