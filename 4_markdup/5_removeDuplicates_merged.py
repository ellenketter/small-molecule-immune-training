# module load python

sample=[
  # '5F1','5F2','5F3',
  #'BG1',
  'BG3',
  #'Fen1',
  'Fen2'
  #,'Fen3','Fluni1',
  #'Fluni2','Fluni3','HC1',
  #'HC2','HC3','HQ1',
  #'HQ2','Myr1','Myr2',
  #'Myr3','Nerol1','Nerol2',
  #'Nerol3','PBS1','PBS2',
  #'PBS3','PBSa1','PBSa2',
  #'PBSa3'
]

sortBams = "/project/lbarreiro/USERS/ellen/KnightMolecules/analysis/5_markdup/5_removeDuplicates.sbatch"

import os
for i in range(1,2):
  os.system("sbatch " + sortBams + " " + sample[i])
  
