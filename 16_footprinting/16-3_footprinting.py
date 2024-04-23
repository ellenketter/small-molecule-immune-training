
sample=[
  '5F',
  'BG',
  'Fen',
  'Fluni',
  'HC',
  'HQ',
  'Myr',
  'Nerol',
  'PBS'
]



foot = "/project/lbarreiro/USERS/ellen/KnightMolecules/analysis/16_footprinting/mergedReplicates/16-3_footprinting.sbatch"

import os
for i in range(8,9):
  os.system("sbatch " + foot + " " + sample[i])
