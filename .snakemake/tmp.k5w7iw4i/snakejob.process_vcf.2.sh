#!/bin/sh
# properties = {"type": "single", "rule": "process_vcf", "local": false, "input": ["output/CTL_samples.txt", "output/FNF_samples.txt"], "output": ["output/geno/cQTL_21_filtered.vcf.gz", "output/geno/pbs_geno/pbs_matched.txt", "output/geno/fnf_geno/fnf_matched.txt"], "wildcards": {}, "params": {"samtoolsVer": "1.21", "raw_vcf": "/work/users/s/e/seyoun/CQTL_AI/Genopipe/output/vcf/CQTL_COA9_10_ALL_qc.vcf.gz", "samples_to_exclude": ""}, "log": ["output/logs/process_vcf.out", "output/logs/process_vcf.err"], "threads": 1, "resources": {}, "jobid": 2, "cluster": {"name": "process_vcf,", "partition": "general", "time": 4320, "cpusPerTask": "1", "memPerCpu": "4G", "nodes": 1, "ntasks": 1, "output": "output/logs_slurm/process_vcf..2.out", "error": "output/logs_slurm/process_vcf..2.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake output/geno/cQTL_21_filtered.vcf.gz --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/VCFpreprocess.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp.k5w7iw4i output/CTL_samples.txt output/FNF_samples.txt --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/ATACconfig.yaml -p --allowed-rules process_vcf --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

