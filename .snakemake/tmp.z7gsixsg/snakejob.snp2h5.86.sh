#!/bin/sh
# properties = {"type": "single", "rule": "snp2h5", "local": false, "input": ["output/wasp/pbs_subset_sample_ids.txt", "output/wasp/fnf_subset_sample_ids.txt"], "output": ["output/wasp/pbs_rep/haplotypes.h5", "output/wasp/pbs_rep/snp_index.h5", "output/wasp/pbs_rep/snp_tab.h5", "output/wasp/fnf_rep/haplotypes.h5", "output/wasp/fnf_rep/snp_index.h5", "output/wasp/fnf_rep/snp_tab.h5"], "wildcards": {}, "params": {"wasp_path": "/nas/longleaf/rhel8/apps/wasp/2023-02/WASP", "chromsize": "/proj/phanstiel_lab/References/chromSizes/hg38_chromSizes_filt_chr.txt", "wasp_ver": "2023-02"}, "log": ["output/logs/pbs_rep_snp2h5.err", "output/logs/fnf_rep_snp2h5.err"], "threads": 1, "resources": {}, "jobid": 86, "cluster": {"name": "snp2h5,", "partition": "general", "time": 4320, "cpusPerTask": "1", "memPerCpu": "4G", "nodes": 1, "ntasks": 1, "output": "output/logs_slurm/snp2h5..86.out", "error": "output/logs_slurm/snp2h5..86.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake output/wasp/pbs_rep/snp_tab.h5 --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/SubsetPreprocess.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp.z7gsixsg output/wasp/pbs_subset_sample_ids.txt output/wasp/fnf_subset_sample_ids.txt --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/ATACconfig.yaml -p --allowed-rules snp2h5 --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

