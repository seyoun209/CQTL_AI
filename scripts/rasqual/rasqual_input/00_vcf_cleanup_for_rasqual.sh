#!/bin/bash

ml plink/1.90b3
ml samtools/1.21
ml vcftools/0.1.15

#This steps are all done from the rename.sh

#bcftools reheader --sample $2 $1 -o $3
#tabix -p vcf $3

#bcftools annotate $3 --rename-chrs /work/users/s/e/seyoun/CQTL_AI/scripts/chrnm.txt -Oz -o $4
#tabix -p vcf $4


#Make bed bim fam files
mkdir -p $1/01.bed
plink --vcf $2 \
        --make-bed \
        --double-id \
        --out $1/01.bed/cqtl

# calculate allele frequency
## at least 3 heterozygous or
## either 3 minor homozygous

mkdir -p $1/02.freq_files

plink --bfile $1/01.bed/cqtl\
        --freqx \
        --out $1/02.freq_files/freq

awk 'NR==1 || ($6 >= 3 || $5 >= 3 )' $1/02.freq_files/freq.frqx > $1/02.freq_files/snp_list.txt;
sed -i '1d' $1/02.freq_files/snp_list.txt;

# filter out snp list
mkdir -p $1/03.filter_out_snps
plink --bfile $1/01.bed/cqtl \
        --extract $1/02.freq_files/snp_list.txt \
        --make-bed \
        --recode 12 \
        --output-missing-genotype 0 \
        --transpose \
        --out $1/03.filter_out_snps/recoded

# change back to vcf file
mkdir -p $1/04.final
plink --bfile $1/03.filter_out_snps/recoded \
        --recode vcf \
        --out $1/04.final/temp_filtered

bgzip $1/04.final/temp_filtered.vcf


bcftools query -l  $1/04.final/temp_filtered.vcf.gz | awk -F_ '{print $0"\t"$1"_"$2"_"$3"_"$4"_"$5}' > $1/04.final/sample_nm.txt

bcftools reheader --sample $1/04.final/sample_nm.txt $1/04.final/temp_filtered.vcf.gz -o $1/04.final/finalFiltered.vcf.gz

tabix -p vcf $1/04.final/finalFiltered.vcf.gz
bcftools stats $1/04.final/finalFiltered.vcf.gz >  $1/04.final/finalFiltered_stats.txt

# Creating PCA

mkdir -p $1/05.pca
plink --vcf $1/04.final/finalFiltered.vcf.gz --pca --out $1/05.pca/cqtl --double-id




