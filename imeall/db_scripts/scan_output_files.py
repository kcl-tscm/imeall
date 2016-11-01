import os
import glob 

jobs = ["00147267160/PotBH/00147267160_v6bxv2_tv0.3bxv0.5/00147267160_v6bxv2_tv0.3bxv0.5_d2.1z",
        "00147267160/PotBH/00147267160_v6bxv2_tv0.1bxv0.2/00147267160_v6bxv2_tv0.1bxv0.2_d2.1z",
        "00147267160/PotBH/00147267160_v6bxv2_tv0.1bxv0.5/00147267160_v6bxv2_tv0.1bxv0.5_d2.1z",
        "00147267160/PotBH/00147267160_v6bxv2_tv0.4bxv0.0/00147267160_v6bxv2_tv0.4bxv0.0_d2.1z",
        "00147267160/PotBH/00147267160_v6bxv2_tv0.1bxv0.3/00147267160_v6bxv2_tv0.1bxv0.3_d2.1z"]

scratch = os.getcwd()
for job in jobs:
  os.chdir(job)
  output_files = sorted(glob.glob('fe*.o*'))
  print job
  os.system('tail -5 {0}'.format(output_files[-1]))
  os.chdir(scratch)
