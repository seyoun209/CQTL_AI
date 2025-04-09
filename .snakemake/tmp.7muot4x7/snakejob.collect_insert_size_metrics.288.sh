#!/bin/sh
# properties = {"type": "single", "rule": "collect_insert_size_metrics", "local": false, "input": ["output/filtered/FACT_23-0023_FNF_synovium_1_filtered.bam"], "output": ["output/metrics/FACT_23-0023_FNF_synovium_1_insert_size_metrics.txt", "output/metrics/FACT_23-0023_FNF_synovium_1_insert_size_histogram.pdf"], "wildcards": {"sampleName": "FACT_23-0023_FNF_synovium_1"}, "params": {"picard_version": "2.26.11", "java_version": "17.0.2"}, "log": ["output/logs/insert_size_metrics_FACT_23-0023_FNF_synovium_1.err", "output/logs/insert_size_metrics_FACT_23-0023_FNF_synovium_1.out"], "threads": 4, "resources": {}, "jobid": 288, "cluster": {"name": "collect_insert_size_metrics,sampleName=FACT_23-0023_FNF_synovium_1", "partition": "general", "time": 4320, "cpusPerTask": "4", "memPerCpu": "6G", "nodes": 1, "ntasks": 1, "output": "output/logs_slurm/collect_insert_size_metrics.sampleName=FACT_23-0023_FNF_synovium_1.288.out", "error": "output/logs_slurm/collect_insert_size_metrics.sampleName=FACT_23-0023_FNF_synovium_1.288.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake output/metrics/FACT_23-0023_FNF_synovium_1_insert_size_metrics.txt --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/ATAC_preprocess.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp.7muot4x7 output/filtered/FACT_23-0023_FNF_synovium_1_filtered.bam --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/ATACconfig.yaml -p --allowed-rules collect_insert_size_metrics --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

