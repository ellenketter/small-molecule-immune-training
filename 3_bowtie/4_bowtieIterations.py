# module load python

sample=[
  'AEK-RK-1s-BG1_S24','AEK-RK-1s-BG3_S25','AEK-RK-1s-Fen1_S21',
  'AEK-RK-1s-Fen2_S22','AEK-RK-1s-Fen3_S23','AEK-RK-1s-Fluni1_S4',
  'AEK-RK-1s-Fluni2_S5','AEK-RK-1s-Fluni3_S6','AEK-RK-1s-HC1_S13',
  'AEK-RK-1s-HC2_S14','AEK-RK-1s-HC3_S15','AEK-RK-1s-HQ1_S19',
  'AEK-RK-1s-HQ2_S20','AEK-RK-1s-Myr1_S10','AEK-RK-1s-Myr2_S11',
  'AEK-RK-1s-Myr3_S12','AEK-RK-1s-Nerol1_S16','AEK-RK-1s-Nerol2_S17',
  'AEK-RK-1s-Nerol3_S18','AEK-RK-1s-PBS1_S1','AEK-RK-1s-PBS2_S2',
  'AEK-RK-1s-PBS3_S3','AEK-RK-1s-PBSa1_S26','AEK-RK-1s-PBSa2_S27',
  'AEK-RK-1s-PBSa3_S28'
]


bowtied = "/project/lbarreiro/USERS/ellen/KnightMolecules/analysis/4_bowtie/4_bowtie.sbatch"

import os
for i in range(12,26):
  os.system("sbatch " + bowtied + " " + sample[i])
  
# out of 28
