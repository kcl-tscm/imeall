import os

files_to_remove = ['bcc_h.xyz.idx',
'defect_cell_relaxed.xyz.idx',
'defect_cell_relaxed.xyz',
'relaxed_cell_removed_defect.xyz',
'relaxed_cell_removed_defect.xyz.idx']


for f in files_to_remove:
    os.remove(f)
