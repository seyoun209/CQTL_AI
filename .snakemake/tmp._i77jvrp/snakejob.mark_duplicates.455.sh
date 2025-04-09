#!/bin/sh
# properties = {"type": "single", "rule": "mark_duplicates", "local": false, "input": ["output/align/FACT_23-0013_CTL_synovium_1_sorted.bam"], "output": ["output/dedup/FACT_23-0013_CTL_synovium_1_dedup.bam", "output/dedup/FACT_23-0013_CTL_synovium_1_dup_metrics.txt", "output/dedup/FACT_23-0013_CTL_synovium_1_dedup.bai"], "wildcards": {"sampleName": "FACT_23-0013_CTL_synovium_1"}, "params": {"picard_version": "2.26.11", "java_version": "17.0.2"}, "log": ["output/logs/mark_duplicates_FACT_23-0013_CTL_synovium_1.err", "output/logs/mark_duplicates_FACT_23-0013_CTL_synovium_1.out"], "threads": 1, "resources": {}, "jobid": 455, "cluster": {"name": "mark_duplicates,sampleName=FACT_23-0013_CTL_synovium_1", "partition": "general", "time": 4320, "cpusPerTask": "1", "memPerCpu": "16G", "nodes": 1, "ntasks": 1, "output": "output/logs_slurm/mark_duplicates.sampleName=FACT_23-0013_CTL_synovium_1.455.out", "error": "output/logs_slurm/mark_duplicates.sampleName=FACT_23-0013_CTL_synovium_1.455.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake output/dedup/FACT_23-0013_CTL_synovium_1_dedup.bam --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/ATAC_preprocess.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp._i77jvrp output/align/FACT_23-0013_CTL_synovium_1_sorted.bam --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/ATACconfig.yaml -p --allowed-rules mark_duplicates --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

