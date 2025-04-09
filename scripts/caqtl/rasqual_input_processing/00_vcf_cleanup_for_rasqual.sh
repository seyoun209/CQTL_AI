#!/bin/bash

ml plink/1.90b3
ml samtools/1.21
ml vcftools/0.1.15

#This steps are all done from the rename.sh

#bcftools reheader --sample $2 $1 -o $3
#tabix -p vcf $3

#bcftools annotate $3 --rename-chrs /work/users/s/e/seyoun/CQTL_AI/scripts/chrnm.txt -Oz -o $4
#tabix -p vcf $4

# Input arguments
# $1 = output directory prefix
# $2 = input VCF file

# Make bed bim fam files
mkdir -p $1/01.bed
plink --vcf $2 \
      --make-bed \
      --double-id \
      --out $1/01.bed/cqtl

# Calculate allele frequency
# At least 2 heterozygous or 2 minor homozygous
mkdir -p $1/02.freq_files
plink --bfile $1/01.bed/cqtl \
      --freqx \
      --out $1/02.freq_files/freq

awk 'NR==1 || ($6 >= 2 || $5 >= 2 )' $1/02.freq_files/freq.frqx > $1/02.freq_files/snp_list.txt
sed -i '1d' $1/02.freq_files/snp_list.txt

# Filter out SNPs and keep the list (no need to recode to VCF yet)
mkdir -p $1/03.filter_out_snps
plink --bfile $1/01.bed/cqtl \
      --extract $1/02.freq_files/snp_list.txt \
      --make-bed \
      --out $1/03.filter_out_snps/filtered

# Extract SNP IDs from the filtered .bim file
awk '{print $2}' $1/03.filter_out_snps/filtered.bim > $1/03.filter_out_snps/filtered_snps.txt

# Filter the original VCF using bcftools to preserve INFO and FORMAT fields
mkdir -p $1/04.final
bcftools view -i 'ID=@'$1/03.filter_out_snps/filtered_snps.txt $2 \
      -o $1/04.final/temp_filtered.vcf -O v

# Compress and index the filtered VCF
bgzip $1/04.final/temp_filtered.vcf
mv $1/04.final/temp_filtered.vcf.gz $1/04.final/temp_filtered.vcf.gz
tabix -p vcf $1/04.final/temp_filtered.vcf.gz

# Update sample names if needed
bcftools query -l $1/04.final/temp_filtered.vcf.gz | awk -F_ '{print $0"\t"$1"_"$2"_"$3"_"$4"_"$5}' > $1/04.final/sample_nm.txt
bcftools reheader --samples $1/04.final/sample_nm.txt $1/04.final/temp_filtered.vcf.gz -o $1/04.final/finalFiltered.vcf.gz

# Index the final VCF
tabix -p vcf $1/04.final/finalFiltered.vcf.gz

# Generate stats
bcftools stats $1/04.final/finalFiltered.vcf.gz > $1/04.final/finalFiltered_stats.txt

# Create PCA
mkdir -p $1/05.pca
plink --vcf $1/04.final/finalFiltered.vcf.gz --pca --out $1/05.pca/cqtl --double-id