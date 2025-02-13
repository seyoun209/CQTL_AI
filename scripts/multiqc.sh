module load multiqc/1.25.1


#Ankle for the verifybamID and the fastqc prior -ATAC (rep and Femur include)
multiqc -f -o /work/users/s/e/seyoun/CQTL_AI/output/QC/multiqc_output/ankle_atac /work/users/s/e/seyoun/CQTL_AI/output/QC/


awk -F'\t' '$7 > 0.05 || $12 > 0.05 {print $0}' /work/users/s/e/seyoun/CQTL_AI/output/QC/multiqc_output/ankle_atac/multiqc_data/multiqc_verifybamid.txt

#everything is below 5%
#Ankle for the verifybamID -RNA (Rep and Femur not included)

multiqc -f -o /work/users/s/e/seyoun/CQTL_AI/rna_output/QC/rna_qc /work/users/s/e/seyoun/CQTL_AI/rna_output/QC

awk -F'\t' '$7 > 0.05 || $12 > 0.05 {print $0}' /work/users/s/e/seyoun/CQTL_AI/rna_output/QC/rna_qc/multiqc_data/multiqc_verifybamid.txt

#CQTL_AM7769_CTL_Ankle_1 sample have the chipmix over 0.92438 
