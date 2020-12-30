#!/bin/bash

# Goal: ancestry admixture analysis of Wellderly subjects
# 10000Genomes VCF is previously preprocessed using the script "prepare_1KG.sh"

# Please, before running the script first edit the following paths
KGDIR=/home/mpinsach/Escriptori/1000Genomes/sequencing/
Wdy_VCF=/home/mpinsach/Escriptori/Wellderly/Wdy_196.vcf
REF_GENOME=/home/mpinsach/Escriptori/ucsc.hg19/ucsc.hg19.fasta
OUTDIR=/home/mpinsach/Escriptori/Ancestry/Wdy/
PRUNED=${OUTDIR}Pruned/
GATK=/home/mpinsach/programes_terminal/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar
ID2POPULATION=/home/mpinsach/Escriptori/scripts/id2population_Wdy.pl
MINMAF=0.05


# Split multiallelic sites
echo Splitting Wellderly VCF
bcftools norm --multiallelics -any ${Wdy_VCF} > ${OUTDIR}Wdy_196_split.vcf 


# Remove Indels
echo Removing indels from Wellderly VCF
java -Xmx20g -jar ${GATK} -R ${REF_GENOME} -T SelectVariants --variant ${OUTDIR}Wdy_196_split.vcf -o ${OUTDIR}Wdy_196_split_SNVs.vcf --selectTypeToInclude SNP 


# Remove alternative positions with "*".
awk '{if ($5!="*") print $0}' ${OUTDIR}Wdy_196_split_SNVs.vcf > ${OUTDIR}Wdy_196_split_SNVs_noasterisk.vcf


# Normalize ID field
echo Normalizing Wellderly VCF
bcftools norm --fasta-ref ${REF_GENOME} --check-ref s ${OUTDIR}Wdy_196_split_SNVs_noasterisk.vcf | awk 'BEGIN{OFS=FS="\t"}{$4=toupper($4);$5=toupper($5); print $0}' | bcftools annotate -x ID --set-id +'%CHROM\_%POS\_%REF\_%ALT' > ${OUTDIR}Wdy_196_split_SNVs_noasterisk_norm.vcf


# Add the population of origin to sample name (ex: Wdy_Samplename)
# This will be used to colour samples in the tSNE plot
perl ${ID2POPULATION} > ${OUTDIR}Wdy_196_split_SNVs_noasterisk_norm_pop.vcf 


# Generate PLINK friently formats
# Here we start working with Wellderly and 1000Genomes VCFs in parallel
plink --vcf ${OUTDIR}Wdy_196_split_SNVs_noasterisk_norm_pop.vcf --keep-allele-order --make-bed --maf ${MINMAF} --out ${OUTDIR}Wdy_196_plink

for i in 3 7 10 11 12
do
plink --vcf ${KGDIR}1000Genomes_chr${i}_split_snvs_target_norm_pop.vcf --keep-allele-order --make-bed --maf ${MINMAF} --out ${OUTDIR}1000Genomes_chr${i}_plink
done 


# LD prune
plink --bfile ${OUTDIR}Wdy_196_plink --indep-pairwise 50 10 0.2 --out ${OUTDIR}Wdy_196_plink.indep.50.10.0.2.g
plink --bfile ${OUTDIR}Wdy_196_plink --extract ${OUTDIR}Wdy_196_plink.indep.50.10.0.2.g.prune.in --make-bed --out ${PRUNED}Wdy_196_pruned

for i in 3 7 10 11 12
do
plink --bfile ${OUTDIR}1000Genomes_chr${i}_plink --indep-pairwise 50 10 0.2 --out ${OUTDIR}1000Genomes_chr${i}_plink.indep.50.10.0.2.g
plink --bfile ${OUTDIR}1000Genomes_chr${i}_plink --extract ${OUTDIR}1000Genomes_chr${i}_plink.indep.50.10.0.2.g.prune.in --make-bed --out ${PRUNED}1000Genomes_chr${i}_pruned
done


# Merge all files
find /home/mpinsach/Escriptori/Ancestry/Wdy/Pruned/ . -name '*pruned.bim' > ${PRUNED}ForMerge.list #Here we create a list of all files ending with .bim
sed -i 's/.bim//g' ${PRUNED}ForMerge.list #Here we replace .bim for .g

plink --merge-list ${PRUNED}ForMerge.list --out ${PRUNED}Merge


# Run PCA 
plink --bfile ${PRUNED}Merge --pca --out ${PRUNED}Wdy_1000Genomes     # PLINK outputs an eigenvec and eigenvalue files.
                                                                      # For ancestry admixture of Wellderly subjects using t-SNE, we use the eigenvec file.
                                                                      # The t-SNE analysis and plot are described in the R script "tSNE_Wdy.R".


