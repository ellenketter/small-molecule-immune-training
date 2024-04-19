# module load python

sample=[
  '5F',
  'BG',
  'Fen',
  'Fluni',
  'HC',
  'HQ',
  'Myr',
  'Nerol'
]

macs2 = "/project/lbarreiro/USERS/ellen/KnightMolecules/analysis/16_footprinting/mergedReplicates/16-6_diffFootprinting_2D.sbatch"

import os
for i in range(0,7):
  os.system("sbatch " + macs2 + " " + sample[i])
  
