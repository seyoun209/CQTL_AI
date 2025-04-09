#!/bin/sh
# properties = {"type": "single", "rule": "signal", "local": false, "input": ["/work/users/s/e/seyoun/CQTL_AI/rna_output/align/FACT_23-0062_FNF_synovium_1_RG.bam"], "output": ["rna_output/signals/synovium/FACT_23-0062_FNF_synovium_1.bw", "rna_output/signals/synovium/norm/FACT_23-0062_FNF_synovium_1.bw"], "wildcards": {"tissue": "synovium", "sampleName": "FACT_23-0062_FNF_synovium_1"}, "params": {"deeptools_ver": "3.5.4", "effective_genomeSize": "2862010428", "NormOption": "CPM"}, "log": ["rna_output/logs/signal_synovium_FACT_23-0062_FNF_synovium_1.err", "rna_output/logs/Normsignal_synovium_FACT_23-0062_FNF_synovium_1.err"], "threads": 1, "resources": {}, "jobid": 38, "cluster": {"name": "signal,sampleName=FACT_23-0062_FNF_synovium_1,tissue=synovium", "partition": "general", "time": 4320, "cpusPerTask": "1", "memPerCpu": "4G", "nodes": 1, "ntasks": 1, "output": "output/logs_slurm/signal.sampleName=FACT_23-0062_FNF_synovium_1,tissue=synovium.38.out", "error": "output/logs_slurm/signal.sampleName=FACT_23-0062_FNF_synovium_1,tissue=synovium.38.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake rna_output/signals/synovium/FACT_23-0062_FNF_synovium_1.bw --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/RNA_signal_track.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp.vhcty_f8 /work/users/s/e/seyoun/CQTL_AI/rna_output/align/FACT_23-0062_FNF_synovium_1_RG.bam --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/RNAconfig.yaml -p --allowed-rules signal --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

