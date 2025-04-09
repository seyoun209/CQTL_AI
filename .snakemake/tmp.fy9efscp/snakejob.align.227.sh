#!/bin/sh
# properties = {"type": "single", "rule": "align", "local": false, "input": ["output/trim/FACT_23-0022_CTL_synovium_1_R1_val_1.fq.gz", "output/trim/FACT_23-0022_CTL_synovium_1_R2_val_2.fq.gz"], "output": ["output/align/FACT_23-0022_CTL_synovium_1_sorted.bam", "output/align/FACT_23-0022_CTL_synovium_1_stats.txt"], "wildcards": {"sampleName": "FACT_23-0022_CTL_synovium_1"}, "params": {"bwa_index": "/users/s/e/seyoun/Ref/genome/bwa_idx/hg38", "bwa_version": "0.7.17", "samtools_version": "1.21"}, "log": ["output/logs/align_FACT_23-0022_CTL_synovium_1.err", "output/logs/align_FACT_23-0022_CTL_synovium_1.out"], "threads": 8, "resources": {}, "jobid": 227, "cluster": {"name": "align,sampleName=FACT_23-0022_CTL_synovium_1", "partition": "general", "time": 4320, "cpusPerTask": "8", "memPerCpu": "12G", "nodes": 1, "ntasks": 1, "output": "output/logs_slurm/align.sampleName=FACT_23-0022_CTL_synovium_1.227.out", "error": "output/logs_slurm/align.sampleName=FACT_23-0022_CTL_synovium_1.227.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake output/align/FACT_23-0022_CTL_synovium_1_sorted.bam --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/ATAC_preprocess.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp.fy9efscp output/trim/FACT_23-0022_CTL_synovium_1_R1_val_1.fq.gz output/trim/FACT_23-0022_CTL_synovium_1_R2_val_2.fq.gz --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/ATACconfig.yaml -p --allowed-rules align --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

