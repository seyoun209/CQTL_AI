#!/bin/sh
# properties = {"type": "single", "rule": "add_read_groups", "local": false, "input": ["output/filtered/CQTL_AM7765_FNF_Ankle_1_filtered.bam", "output/filtered/CQTL_AM7765_FNF_Ankle_1_filtered.bam.bai"], "output": ["output/align/CQTL_AM7765_FNF_Ankle_1_RG.bam", "output/align/CQTL_AM7765_FNF_Ankle_1_RG.bam.bai"], "wildcards": {"sampleName": "CQTL_AM7765_FNF_Ankle_1"}, "params": {"picard_version": "2.26.11", "samtools_version": "1.21"}, "log": ["output/logs/add_read_groups_CQTL_AM7765_FNF_Ankle_1.err", "output/logs/add_read_groups_CQTL_AM7765_FNF_Ankle_1.out"], "threads": 4, "resources": {}, "jobid": 436, "cluster": {"name": "add_read_groups,sampleName=CQTL_AM7765_FNF_Ankle_1", "partition": "general", "time": 4320, "cpusPerTask": "4", "memPerCpu": "6G", "nodes": 1, "ntasks": 1, "output": "output/logs_slurm/add_read_groups.sampleName=CQTL_AM7765_FNF_Ankle_1.436.out", "error": "output/logs_slurm/add_read_groups.sampleName=CQTL_AM7765_FNF_Ankle_1.436.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake output/align/CQTL_AM7765_FNF_Ankle_1_RG.bam --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/ATAC_preprocess.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp.7muot4x7 output/filtered/CQTL_AM7765_FNF_Ankle_1_filtered.bam output/filtered/CQTL_AM7765_FNF_Ankle_1_filtered.bam.bai --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/ATACconfig.yaml -p --allowed-rules add_read_groups --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

