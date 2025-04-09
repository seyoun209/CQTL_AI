#!/bin/bash
  
# Define the list of conditions (adjust if needed)
conditions=(pbs fnf)
windo="100kb"
# Define the list of PCs (pc0 through pc20)
#pcs=(pc0 pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8 pc9 pc10 pc11 pc12 pc13 pc14 pc15 pc16 pc17 pc18 pc19 pc20)
pcs=(pc0 pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8 pc9 pc10 pc11 pc12 pc13)


for cond in "${conditions[@]}"; do
  echo "Processing condition: $cond"
  for pc in "${pcs[@]}"; do
    for chr in {1..22}; do
      echo "  Submitting job for condition: $cond, chromosome: $chr, PC: $pc"
      sbatch /work/users/s/e/seyoun/CQTL_AI/scripts/caqtl/rasqual_output_processing/02_eigenMT_inputPrep.sbatch ${cond} ${pc} ${chr} ${windo}
    done
  done
done
