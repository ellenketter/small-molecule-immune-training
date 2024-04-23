inputs = [
"5f",
"bg",
"fen",
"fluni",
"hc",
"hq",
"myr",
"nerol"
]

homer_annotate = "/project/lbarreiro/USERS/ellen/KnightMolecules/analysis/13_mash/13_homerAnnotatePeaks.sbatch"

import os
for i in range(7,8):
  os.system("sbatch " + homer_annotate + " " + inputs[i])

