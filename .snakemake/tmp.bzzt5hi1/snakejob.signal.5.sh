#!/bin/sh
# properties = {"type": "single", "rule": "signal", "local": false, "input": ["/work/users/s/e/seyoun/CQTL_AI/output/align/FACT_23-0022_CTL_synovium_1_RG.bam"], "output": ["output/signals/synovium/FACT_23-0022_CTL_synovium_1.bw"], "wildcards": {"tissues": "synovium", "sampleName": "FACT_23-0022_CTL_synovium_1"}, "params": {"deeptools_ver": "3.5.4", "binSize": "10", "effective_genomeSize": "2862010428", "NormOption": "RPGC"}, "log": ["output/logs/signal_synovium_FACT_23-0022_CTL_synovium_1.err"], "threads": 1, "resources": {}, "jobid": 5, "cluster": {"name": "signal,sampleName=FACT_23-0022_CTL_synovium_1,tissues=synovium", "partition": "general", "time": 4320, "cpusPerTask": "1", "memPerCpu": "4G", "nodes": 1, "ntasks": 1, "output": "output/logs_slurm/signal.sampleName=FACT_23-0022_CTL_synovium_1,tissues=synovium.5.out", "error": "output/logs_slurm/signal.sampleName=FACT_23-0022_CTL_synovium_1,tissues=synovium.5.err"}}
cd /work/users/s/e/seyoun/CQTL_AI && \
/work/users/s/e/seyoun/CQTL_AI/env/bin/python3 \
-m snakemake output/signals/synovium/FACT_23-0022_CTL_synovium_1.bw --snakefile /work/users/s/e/seyoun/CQTL_AI/snakefiles/signal_track.smk \
--force -j --keep-target-files --keep-remote \
--wait-for-files /work/users/s/e/seyoun/CQTL_AI/.snakemake/tmp.bzzt5hi1 /work/users/s/e/seyoun/CQTL_AI/output/align/FACT_23-0022_CTL_synovium_1_RG.bam --latency-wait 500 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /work/users/s/e/seyoun/CQTL_AI/config/ATACconfig.yaml -p --allowed-rules signal --nocolor --notemp --no-hooks --nolock \
--mode 2  && exit 0 || exit 1

