#!/bin/bash

# Goal: ancestry admixture analysis of GTEx subjects
# 1000 Genomes files used in this script are available at ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/hd_genotype_chip/

# Please, before running the script first edit the following paths
GTEXDIR=/home/mpinsach/Escriptori/GTEx/
KGDIR=/home/mpinsach/Escriptori/1000Genomes/OMNI2.5/
GTEx_VCF=${GTEXDIR}Genotypes/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz
KG_VCF=${KGDIR}/ALL.chip.omni_broad_sanger_combined.20140818.snps.genotypes.vcf.gz
REF_GENOME=/home/mpinsach/Escriptori/human_genome_b37/human_g1k_v37.fasta
TEMPDIR=${GTEXDIR}genotypePCA/tmp
OUTDIR=/home/mpinsach/Escriptori/Ancestry/GTEx/

PREFIX=${OUTDIR}GTEx_wgs_652
PREFIX2=${OUTDIR}1000Genomes_omni
ALLPREFIX=${OUTDIR}GTEx_1000Genomes
MINMAF=0.05

# Remove GL chromosomes from 1000Genomes array
#echo Removing GL chromosomes from 1000 Genomes VCF
#zcat ${KG_VCF} | sed '/^GL/d' | bgzip -c > ${PREFIX2}.noGL.vcf.gz

# Normalize ID field
#echo Normalizing GTEx VCF
#bcftools norm --fasta-ref ${REF_GENOME} --check-ref s ${GTEx_VCF} | awk 'BEGIN{OFS=FS="\t"}{$4=toupper($4);$5=toupper($5);print $0}' | bcftools annotate -x ID --set-id +'%CHROM\_%POS\_%REF\_%ALT' | bgzip -c > ${PREFIX}.norm.vcf.gz
#echo Normalizing 1000Genomes VCF
#bcftools norm --fasta-ref ${REF_GENOME} --check-ref s ${PREFIX2}.noGL.vcf.gz | awk 'BEGIN{OFS=FS="\t"}{$4=toupper($4);$5=toupper($5);print $0}' | bcftools annotate -x ID --set-id +'%CHROM\_%POS\_%REF\_%ALT' | bgzip -c > ${PREFIX2}.norm.vcf.gz


# Convert VCFs to ped format
echo Converting GTEx VCF to ped format
vcftools --temp ${TEMPDIR} --gzvcf ${PREFIX}.norm.vcf.gz --plink --out ${PREFIX} --remove-indels --remove-filtered-all --maf ${MINMAF}
echo Converting 1000Genomes VCF to ped format
vcftools --temp ${TEMPDIR} --gzvcf ${PREFIX2}.norm.vcf.gz --plink --out ${PREFIX2} --remove-indels --remove-filtered-all --maf ${MINMAF}

# Get only biallelic SNPs
echo Getting biallelic SNPs only from GTEx and 1000Genomes
plink --file ${PREFIX} --biallelic-only --out ${PREFIX}.biallelic --make-bed
plink --file ${PREFIX2} --biallelic-only --out ${PREFIX2}.biallelic --make-bed


# Merge files
plink --bfile ${PREFIX}.biallelic --bmerge ${PREFIX2}.biallelic --make-bed --out ${ALLPREFIX}

# LD prune
plink --bfile ${ALLPREFIX} --indep 50 10 2 --out ${ALLPREFIX}.indep.50.10.2.g
plink --bfile ${ALLPREFIX} --exclude ${ALLPREFIX}.indep.50.10.2.g.prune.out --maf ${MINMAF} --out ${ALLPREFIX}.pruned --recode --geno 0.05 --make-bed


# Run PCA 
plink --bfile ${ALLPREFIX}.pruned --pca --out ${ALLPREFIX}.pruned # PLINK outputs an eigenvec and eigenvalue files.
                                                                  # We use the GTEx_1000Genomes.pruned.eigenvec results to plot ancestry admixture of GTEx subjects using ggplot2.
                                        

# Remove intermediate files
rm ${PREFIX}*
rm ${PREFIX2}*
