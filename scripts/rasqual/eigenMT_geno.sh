#/bin/bash

ml plink/2.00a-20220129
ml samtools/1.21

mkdir -p /work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/vcf/recode_012/pbs
mkdir -p /work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/vcf/recode_012/fnf

recode_dir="/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/vcf/recode_012"
base_dir='/work/users/s/e/seyoun/CQTL_AI/output/rasqual_qtl/vcf'


# Process each condition
for CONDITION in fnf pbs; 
do
    echo "Processing condition: ${CONDITION}"
    
    # Create a list of VCF files for this condition
#    VCF_LIST=${base_dir}/${CONDITION}_vcf_list.txt
#    find ${base_dir}/${CONDITION} -name "chr*.${CONDITION}.ASCounts.vcf.gz" |sort -V -k1,1 > ${VCF_LIST}

     COMBINED_VCF=${base_dir}/${CONDITION}/${CONDITION}_combined.vcf



#bcftools concat -f ${VCF_LIST} -o ${COMBINED_VCF} 
#bgzip ${COMBINED_VCF}
#tabix -p vcf ${base_dir}/${CONDITION}/${CONDITION}_combined.vcf.gz

#covert to recode

plink2 --vcf ${base_dir}/${CONDITION}/${CONDITION}_combined.vcf.gz \
	--double-id \
	--recode A-transpose \
	--keep-allele-order \
	--out ${recode_dir}/${CONDITION}/recodeA_${CONDITION}

done




