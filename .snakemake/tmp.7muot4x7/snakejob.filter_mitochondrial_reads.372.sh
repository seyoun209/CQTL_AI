#!/bin/sh
# properties = {"type": "single", "rule": "filter_mitochondrial_reads", "local": false, "input": ["output/dedup/CQTL_AM7755_CTL_Ankle_replicate_dedup.bam", "output/dedup/CQTL_AM7755_CTL_Ankle_replicate_dedup.bam.bai"], "output": ["output/filtered/CQTL_AM7755_CTL_Ankle_replicate_filtered.bam", "output/filtered/CQTL_AM7755_CTL_Ankle_replicate_filtered.bam.bai"], "wildcards": {"sampleName": "CQTL_AM7755_CTL_Ankle_replicate"}, "params": {"samtools_version": "1.21"}, "log": ["output/logs/filter_mitochondrial_CQTL_AM7755_CTL_Ankle_replicate.err", "output/logs/filter_mitochondrial_CQTL_AM7755_CTL_Ankle_replicate.out"], "threads": 4, "resources": {}, "jobid": 372, "cluster": {"name": "filter_mitochondrial_reads,sampleName=CQTL_AM7755_CTL_Ankle_replicate", "partition": "general", "time": 4320, "cpusPerTask": "4", "memPerCpu": "4G", "nodes": 1, "ntasks": 1, "output": "output/logs_slurm/filter_mitochondrial_reads.sampleName=CQTL_AM7755_CTL_Ankle_replicate.372.out", "error": "output/logs_slurm/filter_mitochondrial_reads.sampleName=CQTL_AM7755_CTL_Ankle_replicate.372.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake output/filtered/CQTL_AM7755_CTL_Ankle_replicate_filtered.bam --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/ATAC_preprocess.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp.7muot4x7 output/dedup/CQTL_AM7755_CTL_Ankle_replicate_dedup.bam output/dedup/CQTL_AM7755_CTL_Ankle_replicate_dedup.bam.bai --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/ATACconfig.yaml -p --allowed-rules filter_mitochondrial_reads --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

