#!/bin/bash

# Define the list of conditions (adjust if needed)
conditions=(pbs fnf)

# Define the list of PCs (pc0 through pc20)
#pcs=(pc0 pc1 pc2 pc3 pc4 pc5 pc6 pc7 pc8 pc9 pc10 pc11 pc12 pc13 pc14 pc15 pc16 pc17 pc18 pc19 pc20)
pcs=(pc15 pc16 pc17 pc18 pc19 pc20)

for cond in "${conditions[@]}"; do
  echo "Processing condition: $cond"
  for pc in "${pcs[@]}"; do
    echo "  Submitting job for PC: $pc"
    sbatch /work/users/s/e/seyoun/CQTL_AI/scripts/caqtl/rasqual_output_processing/00_rasqual_combine_output.sbatch ${cond} ${pc}
  done
done

