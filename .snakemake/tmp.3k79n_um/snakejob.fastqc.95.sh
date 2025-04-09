#!/bin/sh
# properties = {"type": "single", "rule": "fastqc", "local": false, "input": ["output/fastq/CQTL_AM7778FE_CTL_Femur_1_R1.fastq.gz", "output/fastq/CQTL_AM7778FE_CTL_Femur_1_R2.fastq.gz"], "output": ["output/QC/CQTL_AM7778FE_CTL_Femur_1_R1_fastqc.zip", "output/QC/CQTL_AM7778FE_CTL_Femur_1_R2_fastqc.zip", "output/QC/CQTL_AM7778FE_CTL_Femur_1_R1_fastqc.html", "output/QC/CQTL_AM7778FE_CTL_Femur_1_R2_fastqc.html"], "wildcards": {"sampleName": "CQTL_AM7778FE_CTL_Femur_1"}, "params": {"dir": "output/QC", "version": "0.12.1"}, "log": ["output/logs/fastqc_CQTL_AM7778FE_CTL_Femur_1.err", "output/logs/fastqc_CQTL_AM7778FE_CTL_Femur_1.out"], "threads": 1, "resources": {}, "jobid": 95, "cluster": {"name": "fastqc,sampleName=CQTL_AM7778FE_CTL_Femur_1", "partition": "general", "time": 4320, "cpusPerTask": "1", "memPerCpu": "4G", "nodes": 1, "ntasks": 1, "output": "output/logs_slurm/fastqc.sampleName=CQTL_AM7778FE_CTL_Femur_1.95.out", "error": "output/logs_slurm/fastqc.sampleName=CQTL_AM7778FE_CTL_Femur_1.95.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake output/QC/CQTL_AM7778FE_CTL_Femur_1_R1_fastqc.zip --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/ATAC_preprocess.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp.3k79n_um output/fastq/CQTL_AM7778FE_CTL_Femur_1_R1.fastq.gz output/fastq/CQTL_AM7778FE_CTL_Femur_1_R2.fastq.gz --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/ATACconfig.yaml -p --allowed-rules fastqc --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

