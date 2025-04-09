#!/bin/sh
# properties = {"type": "single", "rule": "reheader_imputed", "local": false, "input": ["/work/users/s/e/seyoun/CQTL_AI/Genopipe/output/imputation/michigan_imputation_result/chr18.dose.vcf.gz", "output/geno/pbs_geno_other/mappings/sample_name_mapping.txt", "output/geno/fnf_geno_other/mappings/sample_name_mapping.txt"], "output": ["output/geno/pbs_geno_other/wasp_vcf/chr18.pbs.rename.vcf.gz", "output/geno/fnf_geno_other/wasp_vcf/chr18.fnf.rename.vcf.gz"], "wildcards": {"chrom": "18"}, "params": {"chromosomes": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22], "samtools_version": "1.21"}, "log": ["output/logs/imputed_PBS_OTHER_chr18_rename.log", "output/logs/imputed_FNF_OTHER_chr18_rename.log"], "threads": 1, "resources": {}, "jobid": 25, "cluster": {"name": "reheader_imputed,chrom=18", "partition": "general", "time": 4320, "cpusPerTask": "1", "memPerCpu": "4G", "nodes": 1, "ntasks": 1, "output": "output/logs_slurm/reheader_imputed.chrom=18.25.out", "error": "output/logs_slurm/reheader_imputed.chrom=18.25.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake output/geno/pbs_geno_other/wasp_vcf/chr18.pbs.rename.vcf.gz --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/SubsetPreprocess.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp.z7gsixsg /work/users/s/e/seyoun/CQTL_AI/Genopipe/output/imputation/michigan_imputation_result/chr18.dose.vcf.gz output/geno/pbs_geno_other/mappings/sample_name_mapping.txt output/geno/fnf_geno_other/mappings/sample_name_mapping.txt --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/ATACconfig.yaml -p --allowed-rules reheader_imputed --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

