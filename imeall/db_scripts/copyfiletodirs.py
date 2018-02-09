import os
import sys
import glob
import shutil

pattern     = sys.argv[1]
thingtocopy = sys.argv[2]
jobs        = glob.glob(pattern)

print 'file to copy', thingtocopy, 'pattern', pattern
print jobs
for job in jobs:
    shutil.copy(thingtocopy, job)
    print 'Copied', thingtocopy, 'to', job
