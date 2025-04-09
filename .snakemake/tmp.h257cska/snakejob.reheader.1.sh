#!/bin/sh
# properties = {"type": "single", "rule": "reheader", "local": false, "input": ["Genopipe/output/vcf/cQTL_21_filtered.vcf.gz", "Genopipe/output/geno/pbs_geno/pbs_matched.txt", "Genopipe/output/geno/fnf_geno/fnf_matched.txt"], "output": ["Genopipe/output/geno/pbs_geno/pbs.rename_wCHR.vcf.gz", "Genopipe/output/geno/fnf_geno/fnf.rename_wCHR.vcf.gz", "Genopipe/output/geno/pbs_geno/pbs.rename.vcf.gz", "Genopipe/output/geno/fnf_geno/fnf.rename.vcf.gz"], "wildcards": {}, "params": {"samtoolsVer": "1.21"}, "log": ["output/logs/PBS_rename.out", "output/logs/PBS_rename.err", "output/logs/FNF_rename.out", "output/logs/FNF_rename.err"], "threads": 1, "resources": {}, "jobid": 1, "cluster": {"name": "reheader,", "partition": "general", "time": 4320, "cpusPerTask": "1", "memPerCpu": "4G", "nodes": 1, "ntasks": 1, "output": "output/logs_slurm/reheader..1.out", "error": "output/logs_slurm/reheader..1.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake Genopipe/output/geno/pbs_geno/pbs.rename_wCHR.vcf.gz --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/VCFpreprocess.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp.h257cska Genopipe/output/vcf/cQTL_21_filtered.vcf.gz Genopipe/output/geno/pbs_geno/pbs_matched.txt Genopipe/output/geno/fnf_geno/fnf_matched.txt --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/ATACconfig.yaml -p --allowed-rules reheader --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

