# module load python

sample=[
  'AEK-HRK-1s-5F1_S7','AEK-HRK-1s-5F2_S8','AEK-HRK-1s-5F3_S9',
  'AEK-HRK-1s-BG1_S24','AEK-HRK-1s-BG3_S25','AEK-HRK-1s-Fen1_S21',
  'AEK-HRK-1s-Fen2_S22','AEK-HRK-1s-Fen3_S23','AEK-HRK-1s-Fluni1_S4',
  'AEK-HRK-1s-Fluni2_S5','AEK-HRK-1s-Fluni3_S6','AEK-HRK-1s-HC1_S13',
  'AEK-HRK-1s-HC2_S14','AEK-HRK-1s-HC3_S15','AEK-HRK-1s-HQ1_S19',
  'AEK-HRK-1s-HQ2_S20','AEK-HRK-1s-Myr1_S10','AEK-HRK-1s-Myr2_S11',
  'AEK-HRK-1s-Myr3_S12','AEK-HRK-1s-Nerol1_S16','AEK-HRK-1s-Nerol2_S17',
  'AEK-HRK-1s-Nerol3_S18','AEK-HRK-1s-PBS1_S1','AEK-HRK-1s-PBS2_S2',
  'AEK-HRK-1s-PBS3_S3','AEK-HRK-1s-PBSa1_S26','AEK-HRK-1s-PBSa2_S27',
  'AEK-HRK-1s-PBSa3_S28'
]


adapterRemovalNG = "/project/lbarreiro/USERS/ellen/KnightMolecules/analysis/3_NGmerge/3_NGmerge.sbatch"

import os
for i in range(2,28):
  os.system("sbatch " + adapterRemovalNG + " " + sample[i])
