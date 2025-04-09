#!/bin/sh
# properties = {"type": "single", "rule": "signal", "local": false, "input": ["/work/users/s/e/seyoun/CQTL_AI/rna_output/align/CQTL_AM7763_FNF_Ankle_1_RG.bam"], "output": ["rna_output/signals/Ankle/CQTL_AM7763_FNF_Ankle_1.bw", "rna_output/signals/Ankle/norm/CQTL_AM7763_FNF_Ankle_1.bw"], "wildcards": {"tissue": "Ankle", "sampleName": "CQTL_AM7763_FNF_Ankle_1"}, "params": {"deeptools_ver": "3.5.4", "effective_genomeSize": "2862010428", "NormOption": "CPM"}, "log": ["rna_output/logs/signal_Ankle_CQTL_AM7763_FNF_Ankle_1.err", "rna_output/logs/Normsignal_Ankle_CQTL_AM7763_FNF_Ankle_1.err"], "threads": 1, "resources": {}, "jobid": 12, "cluster": {"name": "signal,sampleName=CQTL_AM7763_FNF_Ankle_1,tissue=Ankle", "partition": "general", "time": 4320, "cpusPerTask": "1", "memPerCpu": "4G", "nodes": 1, "ntasks": 1, "output": "output/logs_slurm/signal.sampleName=CQTL_AM7763_FNF_Ankle_1,tissue=Ankle.12.out", "error": "output/logs_slurm/signal.sampleName=CQTL_AM7763_FNF_Ankle_1,tissue=Ankle.12.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake rna_output/signals/Ankle/CQTL_AM7763_FNF_Ankle_1.bw --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/RNA_signal_track.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp.vhcty_f8 /work/users/s/e/seyoun/CQTL_AI/rna_output/align/CQTL_AM7763_FNF_Ankle_1_RG.bam --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/RNAconfig.yaml -p --allowed-rules signal --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

