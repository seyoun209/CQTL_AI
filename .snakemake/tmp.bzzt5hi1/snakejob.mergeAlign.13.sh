#!/bin/sh
# properties = {"type": "single", "rule": "mergeAlign", "local": false, "input": ["/work/users/s/e/seyoun/CQTL_AI/output/align/FACT_23-0013_CTL_synovium_1_RG.bam", "/work/users/s/e/seyoun/CQTL_AI/output/align/FACT_23-0016_CTL_synovium_1_RG.bam", "/work/users/s/e/seyoun/CQTL_AI/output/align/FACT_23-0022_CTL_synovium_1_RG.bam", "/work/users/s/e/seyoun/CQTL_AI/output/align/FACT_23-0023_CTL_synovium_1_RG.bam", "/work/users/s/e/seyoun/CQTL_AI/output/align/FACT_23-0025_CTL_synovium_1_RG.bam"], "output": ["output/siganls/synovium/mergedAlign/CTL_sorted.bam", "output/siganls/synovium/mergedAlign/CTL_sorted.bam.bai", "output/siganls/synovium/mergedAlign/CTL_stats.txt"], "wildcards": {"tissues": "synovium", "cond": "CTL"}, "params": {"samtools_version": "1.21"}, "log": ["output/logs/mergedAligned_synovium_CTL.err"], "threads": 1, "resources": {}, "jobid": 13, "cluster": {"name": "mergeAlign,cond=CTL,tissues=synovium", "partition": "general", "time": 4320, "cpusPerTask": "1", "memPerCpu": "4G", "nodes": 1, "ntasks": 1, "output": "output/logs_slurm/mergeAlign.cond=CTL,tissues=synovium.13.out", "error": "output/logs_slurm/mergeAlign.cond=CTL,tissues=synovium.13.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake output/siganls/synovium/mergedAlign/CTL_sorted.bam --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/signal_track.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp.bzzt5hi1 /work/users/s/e/seyoun/CQTL_AI/output/align/FACT_23-0013_CTL_synovium_1_RG.bam /work/users/s/e/seyoun/CQTL_AI/output/align/FACT_23-0016_CTL_synovium_1_RG.bam /work/users/s/e/seyoun/CQTL_AI/output/align/FACT_23-0022_CTL_synovium_1_RG.bam /work/users/s/e/seyoun/CQTL_AI/output/align/FACT_23-0023_CTL_synovium_1_RG.bam /work/users/s/e/seyoun/CQTL_AI/output/align/FACT_23-0025_CTL_synovium_1_RG.bam --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/ATACconfig.yaml -p --allowed-rules mergeAlign --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

