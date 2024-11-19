#!/bin/bash

input_gtf="/users/s/e/seyoun/Ref/genome/gencode.v45.annotation.gtf"
output_bed="/users/s/e/seyoun/Ref/genome/gencode.v45.annotation.tss.bed"

# Extract TSS from the GTF file
awk 'BEGIN {FS="\t"; OFS="\t"} $3 == "transcript" {
    split($9, attributes, ";");
    for (i in attributes) {
        if (attributes[i] ~ /transcript_id/) {
            gsub(/"/, "", attributes[i]);
            gsub(/transcript_id /, "", attributes[i]);
            transcript_id = attributes[i];
        }
        if (attributes[i] ~ /gene_name/) {
            gsub(/"/, "", attributes[i]);
            gsub(/gene_name /, "", attributes[i]);
            gene_name = attributes[i];
        }
    }
    if ($7 == "+") {
        print $1, $4-1, $4, transcript_id, gene_name, $7;
    } else if ($7 == "-") {
        print $1, $5-1, $5, transcript_id, gene_name, $7;
    }
}' $input_gtf | sort -k1,1 -k2,2n > $output_bed

echo "TSS BED file created: $output_bed"

