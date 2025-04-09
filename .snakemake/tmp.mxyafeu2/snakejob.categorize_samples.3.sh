#!/bin/sh
# properties = {"type": "single", "rule": "categorize_samples", "local": false, "input": [], "output": ["output/CTL_samples.txt", "output/FNF_samples.txt"], "wildcards": {}, "params": {}, "log": ["output/logs/categorize_samples.out", "output/logs/categorize_samples.err"], "threads": 1, "resources": {}, "jobid": 3, "cluster": {"name": "categorize_samples,", "partition": "general", "time": 4320, "cpusPerTask": "1", "memPerCpu": "4G", "nodes": 1, "ntasks": 1, "output": "output/logs_slurm/categorize_samples..3.out", "error": "output/logs_slurm/categorize_samples..3.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake output/CTL_samples.txt --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/VCFpreprocess.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp.mxyafeu2 --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/ATACconfig.yaml -p --allowed-rules categorize_samples --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

