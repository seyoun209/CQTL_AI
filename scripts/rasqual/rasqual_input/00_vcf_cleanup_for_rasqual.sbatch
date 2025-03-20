#!/bin/bash
#SBATCH -J vcfcleanup
#SBATCH -p general
#SBATCH -n 2
#SBATCH -N 1
#SBATCH --mem=15g
#SBATCH -t 1-00:00:00
#SBATCH -e ../../output/logs/fnf_vcfcleanup_e.%j
#SBATCH -o ../output/logs/fnf_vcfcleanup_o.%j


#This is PBS
sh ./vcf_cleanup_for_rasqual.sh /work/users/s/e/seyoun/CQTL_AI/output/geno/pbs_geno /work/users/s/e/seyoun/CQTL_AI/output/geno/pbs_geno/pbs.rename_wCHR.vcf.gz 

#This is FNF
sh ./vcf_cleanup_for_rasqual.sh /work/users/s/e/seyoun/CQTL_AI/output/geno/fnf_geno /work/users/s/e/seyoun/CQTL_AI/output/geno/fnf_geno/fnf.rename_wCHR.vcf.gz

