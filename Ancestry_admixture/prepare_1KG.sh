#!/bin/bash

# Goal: pre-process 1000 Genomes VCF files to be used for ancestry admixture analysis of Brugada syndrome and Wellderly subjects (scripts "ancestry_BrS.sh" & "ancestry_Wdy.sh")
# 1000 Genomes files used in this script are available at ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/

# Please, before running the script first edit the following paths
KGDIR=/home/mpinsach/Escriptori/1000Genomes/sequencing
TARGET_REGIONS=/home/mpinsach/Escriptori/target_regions/RegSeq_numbered_nochr.bed
GATK=/home/mpinsach/programes_terminal/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar
REF_GENOME=/home/mpinsach/Escriptori/human_genome_b37/human_g1k_v37.fasta
ID2POPULATION=/home/mpinsach/Escriptori/scripts/id2population_1KG.pl


# Split multiallelic sites
for i in 3 7 10 11 12
do
echo Splitting chr${i}
bcftools norm --multiallelics -any ${KGDIR}/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | bgzip -c > ${KGDIR}/1000Genomes_chr${i}_split.vcf.gz
done

# Remove all variants that are not SNVs
for i in 3 7 10 11 12
do
echo Keeping only SNVs from chr${i}
vcftools --gzvcf ${KGDIR}/1000Genomes_chr${i}_split.vcf.gz --remove-indels --remove-filtered-all --recode --recode-INFO-all --out ${KGDIR}/1000Genomes_chr${i}_split_snvs
done

# Intersect with target cis-regulatory regions
for i in 3 7 10 11 12
do
echo Intersecting chr${i}
bedtools intersect -header -wa -a ${KGDIR}/1000Genomes_chr${i}_split_snvs.recode.vcf -b /home/mpinsach/Escriptori/target_regions/RegSeq_numbered_nochr.bed > ${KGDIR}/1000Genomes_chr${i}_split_snvs_target.vcf
done 

# Normalize ID filed
for i in 3 7 10 11 12
do
echo Normalizing chr${i}
bcftools norm --fasta-ref ${REF_GENOME} --check-ref s ${KGDIR}/1000Genomes_chr${i}_split_snvs_target.vcf | awk 'BEGIN{OFS=FS="\t"}{$4=toupper($4);$5=toupper($5);print $0}' | bcftools annotate -x ID --set-id +'%CHROM\_%POS\_%REF\_%ALT' | awk 'BEGIN {OFS=FS="\t"}{if ($0!~"^#"){$3="chr"$3; print $0} else print $0}' > ${KGDIR}/1000Genomes_chr${i}_split_snvs_target_norm.vcf
done

# Add the population of origin to sample name (ex: AFR_HG01885)
# This will be used to colour samples in the tSNE plot
perl ${ID2POPULATION}
