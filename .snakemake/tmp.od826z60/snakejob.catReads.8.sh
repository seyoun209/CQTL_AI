#!/bin/sh
# properties = {"type": "single", "rule": "catReads", "local": false, "input": ["/proj/phanstiel_lab/Data/raw/CQTL/atac/02132023_nextseq_i7/CQTL_AM7711_A_FNF_18_1_S2_R1_001.fastq.gz", "/proj/phanstiel_lab/Data/raw/CQTL/atac/02132023_nextseq_i7/CQTL_AM7711_A_FNF_18_1_S2_R2_001.fastq.gz"], "output": ["output/fastq/CQTL_AM7711_FNF_Ankle_1_R1.fastq.gz", "output/fastq/CQTL_AM7711_FNF_Ankle_1_R2.fastq.gz"], "wildcards": {"sampleName": "CQTL_AM7711_FNF_Ankle_1"}, "params": {}, "log": ["output/logs/CQTL_AM7711_FNF_Ankle_1_R1_catReads.err", "output/logs/CQTL_AM7711_FNF_Ankle_1_R2_catReads.err"], "threads": 1, "resources": {}, "jobid": 8, "cluster": {"name": "catReads,sampleName=CQTL_AM7711_FNF_Ankle_1", "partition": "general", "time": 4320, "cpusPerTask": "1", "memPerCpu": "4G", "nodes": 1, "ntasks": 1, "output": "output/logs_slurm/catReads.sampleName=CQTL_AM7711_FNF_Ankle_1.8.out", "error": "output/logs_slurm/catReads.sampleName=CQTL_AM7711_FNF_Ankle_1.8.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake output/fastq/CQTL_AM7711_FNF_Ankle_1_R1.fastq.gz --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/ATAC_preprocess.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp.od826z60 /proj/phanstiel_lab/Data/raw/CQTL/atac/02132023_nextseq_i7/CQTL_AM7711_A_FNF_18_1_S2_R1_001.fastq.gz /proj/phanstiel_lab/Data/raw/CQTL/atac/02132023_nextseq_i7/CQTL_AM7711_A_FNF_18_1_S2_R2_001.fastq.gz --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/ATACconfig.yaml -p --allowed-rules catReads --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

