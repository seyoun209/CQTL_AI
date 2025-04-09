#!/bin/sh
# properties = {"type": "single", "rule": "mergeAlign", "local": false, "input": ["/work/users/s/e/seyoun/CQTL_AI/rna_output/align/CQTL_AM7778FE_FNF_Femur_1_RG.bam"], "output": ["rna_output/signals/Femur/mergedAlign/CTL_sorted.bam", "rna_output/signals/Femur/mergedAlign/CTL_sorted.bam.bai", "rna_output/signals/Femur/mergedAlign/CTL_stats.txt"], "wildcards": {"tissue": "Femur", "cond": "CTL"}, "params": {"samtools_version": "1.21"}, "log": ["rna_output/logs/mergedAligned_Femur_CTL.err"], "threads": 1, "resources": {}, "jobid": 63, "cluster": {"name": "mergeAlign,cond=CTL,tissue=Femur", "partition": "general", "time": 4320, "cpusPerTask": "1", "memPerCpu": "8G", "nodes": 1, "ntasks": 1, "output": "output/logs_slurm/mergeAlign.cond=CTL,tissue=Femur.63.out", "error": "output/logs_slurm/mergeAlign.cond=CTL,tissue=Femur.63.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake rna_output/signals/Femur/mergedAlign/CTL_sorted.bam --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/RNA_signal_track.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp.vhcty_f8 /work/users/s/e/seyoun/CQTL_AI/rna_output/align/CQTL_AM7778FE_FNF_Femur_1_RG.bam --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/RNAconfig.yaml -p --allowed-rules mergeAlign --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

