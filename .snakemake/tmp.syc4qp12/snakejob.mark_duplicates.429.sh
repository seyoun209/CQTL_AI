#!/bin/sh
# properties = {"type": "single", "rule": "mark_duplicates", "local": false, "input": ["output/align/CQTL_AM7755_CTL_Ankle_1_sorted.bam"], "output": ["output/dedup/CQTL_AM7755_CTL_Ankle_1_dedup.bam", "output/dedup/CQTL_AM7755_CTL_Ankle_1_dup_metrics.txt", "output/dedup/CQTL_AM7755_CTL_Ankle_1_dedup.bai"], "wildcards": {"sampleName": "CQTL_AM7755_CTL_Ankle_1"}, "params": {"picard_version": "2.26.11", "java_version": "17.0.2"}, "log": ["output/logs/mark_duplicates_CQTL_AM7755_CTL_Ankle_1.err", "output/logs/mark_duplicates_CQTL_AM7755_CTL_Ankle_1.out"], "threads": 1, "resources": {}, "jobid": 429, "cluster": {"name": "mark_duplicates,sampleName=CQTL_AM7755_CTL_Ankle_1", "partition": "general", "time": 4320, "cpusPerTask": "1", "memPerCpu": "16G", "nodes": 1, "ntasks": 4, "output": "output/logs_slurm/mark_duplicates.sampleName=CQTL_AM7755_CTL_Ankle_1.429.out", "error": "output/logs_slurm/mark_duplicates.sampleName=CQTL_AM7755_CTL_Ankle_1.429.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake output/dedup/CQTL_AM7755_CTL_Ankle_1_dedup.bam --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/ATAC_preprocess.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp.syc4qp12 output/align/CQTL_AM7755_CTL_Ankle_1_sorted.bam --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/ATACconfig.yaml -p --allowed-rules mark_duplicates --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

