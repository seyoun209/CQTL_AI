#!/bin/sh
# properties = {"type": "single", "rule": "extract_sample_ids", "local": false, "input": ["output/geno/pbs_geno_other/CTL_samples_with_replacements.txt", "output/geno/fnf_geno_other/FNF_samples_with_replacements.txt"], "output": ["output/wasp/pbs_subset_sample_ids.txt", "output/wasp/fnf_subset_sample_ids.txt"], "wildcards": {}, "params": {}, "log": ["output/logs/extract_sample_ids_pbs_subset.log", "output/logs/extract_sample_ids_fnf_subset.log"], "threads": 1, "resources": {}, "jobid": 87, "cluster": {"name": "extract_sample_ids,", "partition": "general", "time": 4320, "cpusPerTask": "1", "memPerCpu": "4G", "nodes": 1, "ntasks": 1, "output": "output/logs_slurm/extract_sample_ids..87.out", "error": "output/logs_slurm/extract_sample_ids..87.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake output/wasp/pbs_subset_sample_ids.txt --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/SubsetPreprocess.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp.z7gsixsg output/geno/pbs_geno_other/CTL_samples_with_replacements.txt output/geno/fnf_geno_other/FNF_samples_with_replacements.txt --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/ATACconfig.yaml -p --allowed-rules extract_sample_ids --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

