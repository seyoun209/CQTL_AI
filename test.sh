module load java/17.0.2
module load picard/2.26.11
ml r/4.4.0
java -jar /nas/longleaf/apps/picard/2.26.11/picard-2.26.11/picard.jar CollectInsertSizeMetrics \
    I=output/filtered/CQTL_AM7717_FNF_Ankle_1_filtered.bam \
    H=output/metrics/CQTL_AM7717_FNF_Ankle_1_insert_size_histogram.pdf \
    M=0.05 \
    ASSUME_SORTED=true \
    O=output/metrics/CQTL_AM7717_FNF_Ankle_1_insert_size_metrics.txt
