#!/bin/sh
# properties = {"type": "single", "rule": "mergeAlign", "local": false, "input": ["/work/users/s/e/seyoun/CQTL_AI/rna_output/align/CQTL_AM7778_FNF_Ankle_1_RG.bam"], "output": ["rna_output/signals/Femur/mergedAlign/FNF_sorted.bam", "rna_output/signals/Femur/mergedAlign/FNF_sorted.bam.bai", "rna_output/signals/Femur/mergedAlign/FNF_stats.txt"], "wildcards": {"tissue": "Femur", "cond": "FNF"}, "params": {"samtools_version": "1.21"}, "log": ["rna_output/logs/mergedAligned_Femur_FNF.err"], "threads": 1, "resources": {}, "jobid": 64, "cluster": {"name": "mergeAlign,cond=FNF,tissue=Femur", "partition": "general", "time": 4320, "cpusPerTask": "1", "memPerCpu": "8G", "nodes": 1, "ntasks": 1, "output": "output/logs_slurm/mergeAlign.cond=FNF,tissue=Femur.64.out", "error": "output/logs_slurm/mergeAlign.cond=FNF,tissue=Femur.64.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake rna_output/signals/Femur/mergedAlign/FNF_sorted.bam --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/RNA_signal_track.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp.vhcty_f8 /work/users/s/e/seyoun/CQTL_AI/rna_output/align/CQTL_AM7778_FNF_Ankle_1_RG.bam --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/RNAconfig.yaml -p --allowed-rules mergeAlign --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

