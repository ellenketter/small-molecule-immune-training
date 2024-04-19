inputs = [
# "5f",
"bg",
"fen",
"fluni",
"hc",
# "hq",
"myr"
# "nerol",
# "background"
]

homer_annotate = "/project/lbarreiro/USERS/ellen/KnightMolecules/analysis/11_annotatePeaks/11_homerAnnotatePeaks.sbatch"

import os
for i in range(1,5):
  os.system("sbatch " + homer_annotate + " " + inputs[i])

