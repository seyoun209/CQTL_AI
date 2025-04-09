#!/bin/sh
# properties = {"type": "single", "rule": "mergeForwardSignal", "local": false, "input": ["rna_output/signals/synovium/mergedAlign/CTL_sorted.bam"], "output": ["rna_output/signals/synovium/strand/CTL_fwd.bw", "rna_output/signals/synovium/strand/CTL_norm_fwd.bw"], "wildcards": {"tissue": "synovium", "cond": "CTL"}, "params": {"deeptools_ver": "3.5.4", "effective_genomeSize": "2862010428", "NormOption": "CPM"}, "log": ["rna_output/logs/mergeSignal_forward_synovium_CTL.err", "rna_output/logs/Norm_mergeSignal_forward_synovium_CTL.err"], "threads": 1, "resources": {}, "jobid": 57, "cluster": {"name": "mergeForwardSignal,cond=CTL,tissue=synovium", "partition": "general", "time": 4320, "cpusPerTask": "1", "memPerCpu": "4G", "nodes": 1, "ntasks": 1, "output": "output/logs_slurm/mergeForwardSignal.cond=CTL,tissue=synovium.57.out", "error": "output/logs_slurm/mergeForwardSignal.cond=CTL,tissue=synovium.57.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake rna_output/signals/synovium/strand/CTL_fwd.bw --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/RNA_signal_track.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp.vhcty_f8 rna_output/signals/synovium/mergedAlign/CTL_sorted.bam --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/RNAconfig.yaml -p --allowed-rules mergeForwardSignal --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

