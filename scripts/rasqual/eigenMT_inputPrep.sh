#!/bin/bash

# Conditions, chromosomes, and PCs
conditions=("pbs" "fnf")
chromosomes=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 "X")
pcs=(1 2 3 4 5 6 7 8 9 10)

# Base directory for logs
mkdir -p ./logs

# Counter for total jobs
total_jobs=0

# Nested loops to create and submit jobs for each combination
for condition in "${conditions[@]}"; do
    for chr in "${chromosomes[@]}"; do
        for pc in "${pcs[@]}"; do
            # Create job name
            job_name="eigenMT_${condition}_chr${chr}_PC${pc}"
            
            # Output and error log files
            output_log="./logs/${job_name}_%j.out"
            error_log="./logs/${job_name}_%j.err"
            
            # Submit the job
            sbatch <<EOF
#!/bin/bash
#SBATCH -J $job_name
#SBATCH -p general
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --mem=10G
#SBATCH -t 2:00:00
#SBATCH -o $output_log
#SBATCH -e $error_log

# Load R module if needed
module load r/4.4.0

# Run R script with parameters
Rscript eigenMT_inputPrep.R $condition $pc $chr
EOF
            
            # Increment job counter
            ((total_jobs++))
            
            # Optional: Add a small delay to prevent overwhelming the scheduler
            sleep 0.1
        done
    done
done
echo "Submitted $total_jobs jobs"

echo "Submitted $total_jobs jobs"
