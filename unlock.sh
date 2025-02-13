#! /bin/bash -login
  
## Exit if any command fails
set -e

## Load required modules
module load python/3.9.6

## Activate virtual environment
source env/bin/activate

case $1 in

    '-h' | '--help' | ' ' | '')
            echo -e '\e[31mSpecify which workflow to unlock (i.e. run_RNAprocessing)'
            exit 2
            ;;
        'run_ATAC_preprocessing')
            ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s snakefiles/ATAC_preprocess.smk --configfile "config/ATACconfig.yaml"  --cluster-config "config/cluster.yaml" --cluster "sbatch -J {cluster.name} -p {cluster.partition} -t {cluster.time} -c {cluster.cpusPerTask} --mem-per-cpu={cluster.memPerCpu} -N {cluster.nodes} --output {cluster.output} --error {cluster.error} --parsable" --cluster-status snakefiles/utils/status.py
            ;;
    'run_VCFpreprocessing')
            ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s snakefiles/VCFpreprocess.smk --configfile "config/ATACconfig.yaml"  --cluster-config "config/cluster.yaml" --cluster "sbatch -J {cluster.name} -p {cluster.partition} -t {cluster.time} -c {cluster.cpusPerTask} --mem-per-cpu={cluster.memPerCpu} -N {cluster.nodes} --output {cluster.output} --error {cluster.error} --parsable" --cluster-status snakefiles/utils/status.py
            ;;
    'run_SubsetPreprocessing')
            ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s snakefiles/SubsetPreprocess.smk --configfile "config/ATACconfig.yaml"  --cluster-config "config/cluster.yaml" --cluster "sbatch -J {cluster.name} -p {cluster.partition} -t {cluster.time} -c {cluster.cpusPerTask} --mem-per-cpu={cluster.memPerCpu} -N {cluster.nodes} --output {cluster.output} --error {cluster.error} --parsable" --cluster-status snakefiles/utils/status.py
            ;;
	    'rna' | 'run_RNA_preprocessing')
	    ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s snakefiles/RNA_preprocess.smk --configfile "config/RNAconfig.yaml"  --cluster-config "config/cluster.yaml" --cluster "sbatch -J {cluster.name} -p {cluster.partition} -t {cluster.time} -c {cluster.cpusPerTask} --mem-per-cpu={cluster.memPerCpu} -N {cluster.nodes} --output {cluster.output} --error {cluster.error} --parsable" --cluster-status snakefiles/utils/status.py
            ;;
    'wasp' | 'run_WASP_ATAC')
            ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s snakefiles/WASP_atac.smk --configfile "config/ATACconfig.yaml"  --cluster-config "config/cluster.yaml" --cluster "sbatch -J {cluster.name} -p {cluster.partition} -t {cluster.time} -c {cluster.cpusPerTask} --mem-per-cpu={cluster.memPerCpu} -N {cluster.nodes} --output {cluster.output} --error {cluster.error} --parsable" --cluster-status snakefiles/utils/status.py
            ;;
    'signal' | 'run_signalTrack')
            ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s snakefiles/signal_track.smk --configfile "config/ATACconfig.yaml"  --cluster-config "config/cluster.yaml" --cluster "sbatch -J {cluster.name} -p {cluster.partition} -t {cluster.time} -c {cluster.cpusPerTask} --mem-per-cpu={cluster.memPerCpu} -N {cluster.nodes} --output {cluster.output} --error {cluster.error} --parsable" --cluster-status snakefiles/utils/status.py
            ;;
    'peak' | 'run_peakcalling')
            ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s snakefiles/peak_counting.smk --configfile "config/ATACconfig.yaml"  --cluster-config "config/cluster.yaml" --cluster "sbatch -J {cluster.name} -p {cluster.partition} -t {cluster.time} -c {cluster.cpusPerTask} --mem-per-cpu={cluster.memPerCpu} -N {cluster.nodes} --output {cluster.output} --error {cluster.error} --parsable" --cluster-status snakefiles/utils/status.py
            ;;
    'rna_signal' | 'run_RNAsignal')
            ## Unlock snakemake workflow
            snakemake -j 1 --unlock -s snakefiles/RNA_signal_track.smk --configfile "config/RNAconfig.yaml"  --cluster-config "config/cluster.yaml" --cluster "sbatch -J {cluster.name} -p {cluster.partition} -t {cluster.time} -c {cluster.cpusPerTask} --mem-per-cpu={cluster.memPerCpu} -N {cluster.nodes} --output {cluster.output} --error {cluster.error} --parsable" --cluster-status snakefiles/utils/status.py
            ;;
esac

## Deactivate virtual environment
deactivate

## Success message
echo -e "\033[0;32mDirectory unlocked, ready to rerun."
