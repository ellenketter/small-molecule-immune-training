inputs = [
"annotateInput_DOWN_5f.txt",
"annotateInput_DOWN_bg.txt",
"annotateInput_DOWN_fen.txt",
"annotateInput_DOWN_fluni.txt",
"annotateInput_DOWN_hc.txt",
"annotateInput_DOWN_hq.txt",
"annotateInput_DOWN_myr.txt",
"annotateInput_DOWN_nerol.txt",
"annotateInput_UP_5f.txt",
"annotateInput_UP_bg.txt",
"annotateInput_UP_fen.txt",
"annotateInput_UP_fluni.txt",
"annotateInput_UP_hc.txt",
"annotateInput_UP_hq.txt",
"annotateInput_UP_myr.txt",
"annotateInput_UP_nerol.txt"
]


background = [
"annotateInput.txt",
"annotateInput.txt",
"annotateInput.txt",
"annotateInput.txt",
"annotateInput.txt",
"annotateInput.txt",
"annotateInput.txt",
"annotateInput.txt",
"annotateInput.txt",
"annotateInput.txt",
"annotateInput.txt",
"annotateInput.txt",
"annotateInput.txt",
"annotateInput.txt",
"annotateInput.txt",
"annotateInput.txt"
]


output = [

"tfMotifOutput_DOWN_5f",
"tfMotifOutput_DOWN_bg",
"tfMotifOutput_DOWN_fen",
"tfMotifOutput_DOWN_fluni",
"tfMotifOutput_DOWN_hc",
"tfMotifOutput_DOWN_hq",
"tfMotifOutput_DOWN_myr",
"tfMotifOutput_DOWN_nerol",
"tfMotifOutput_UP_5f",
"tfMotifOutput_UP_bg",
"tfMotifOutput_UP_fen",
"tfMotifOutput_UP_fluni",
"tfMotifOutput_UP_hc",
"tfMotifOutput_UP_hq",
"tfMotifOutput_UP_myr",
"tfMotifOutput_UP_nerol"
]



homer_annotate = "/project/lbarreiro/USERS/ellen/KnightMolecules/analysis/15_transcriptionFactorMotifs/15_homer_motif.sbatch"

import os
for i in range(0,16):
  os.system("sbatch " + homer_annotate + " " + inputs[i] + " " + background[i] + " " +  output[i])
