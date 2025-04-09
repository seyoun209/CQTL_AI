#!/bin/sh
# properties = {"type": "single", "rule": "trim", "local": false, "input": ["output/fastq/CQTL_AM7711_FNF_Ankle_1_R1.fastq.gz", "output/fastq/CQTL_AM7711_FNF_Ankle_1_R2.fastq.gz"], "output": ["output/trim/CQTL_AM7711_FNF_Ankle_1_R1_val_1.fq.gz", "output/trim/CQTL_AM7711_FNF_Ankle_1_R2_val_2.fq.gz", "output/trim/CQTL_AM7711_FNF_Ankle_1_R1.fastq.gz_trimming_report.txt", "output/trim/CQTL_AM7711_FNF_Ankle_1_R2.fastq.gz_trimming_report.txt"], "wildcards": {"sampleName": "CQTL_AM7711_FNF_Ankle_1"}, "params": {"version": "0.6.7"}, "log": ["output/logs/trim_CQTL_AM7711_FNF_Ankle_1.err", "output/logs/trim_CQTL_AM7711_FNF_Ankle_1.out"], "threads": 4, "resources": {}, "jobid": 124, "cluster": {"name": "trim,sampleName=CQTL_AM7711_FNF_Ankle_1", "partition": "general", "time": 4320, "cpusPerTask": "4", "memPerCpu": "4G", "nodes": 1, "ntasks": 1, "output": "output/logs_slurm/trim.sampleName=CQTL_AM7711_FNF_Ankle_1.124.out", "error": "output/logs_slurm/trim.sampleName=CQTL_AM7711_FNF_Ankle_1.124.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake output/trim/CQTL_AM7711_FNF_Ankle_1_R1.fastq.gz_trimming_report.txt --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/ATAC_preprocess.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp.ce1meoxy output/fastq/CQTL_AM7711_FNF_Ankle_1_R1.fastq.gz output/fastq/CQTL_AM7711_FNF_Ankle_1_R2.fastq.gz --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/ATACconfig.yaml -p --allowed-rules trim --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

