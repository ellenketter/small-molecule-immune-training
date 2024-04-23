# module load python

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

macs2 = "/project/lbarreiro/USERS/ellen/KnightMolecules/analysis/16_footprinting/mergedReplicates/16-2_peakCalling.sbatch"

import os
for i in range(8,9):
  os.system("sbatch " + macs2 + " " + sample[i])
  
