#!/bin/sh
# properties = {"type": "single", "rule": "mergeReverseSignal", "local": false, "input": ["rna_output/signals/synovium/mergedAlign/FNF_sorted.bam"], "output": ["rna_output/signals/synovium/strand/FNF_rev.bw", "rna_output/signals/synovium/strand/FNF_norm_rev.bw"], "wildcards": {"tissue": "synovium", "cond": "FNF"}, "params": {"deeptools_ver": "3.5.4", "effective_genomeSize": "2862010428", "NormOption": "CPM"}, "log": ["rna_output/logs/mergeSignal_reverse_synovium_FNF.err", "rna_output/logs/Norm_mergeSignal_reverse_synovium_FNF.err"], "threads": 1, "resources": {}, "jobid": 60, "cluster": {"name": "mergeReverseSignal,cond=FNF,tissue=synovium", "partition": "general", "time": "11-00:00:00", "cpusPerTask": "1", "memPerCpu": "16G", "nodes": 1, "ntasks": 1, "output": "output/logs_slurm/mergeReverseSignal.cond=FNF,tissue=synovium.60.out", "error": "output/logs_slurm/mergeReverseSignal.cond=FNF,tissue=synovium.60.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake rna_output/signals/synovium/strand/FNF_rev.bw --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/RNA_signal_track.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp.vhcty_f8 rna_output/signals/synovium/mergedAlign/FNF_sorted.bam --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/RNAconfig.yaml -p --allowed-rules mergeReverseSignal --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

