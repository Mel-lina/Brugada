#!/usr/bin/perl

# Goal: Variant calling of Brugada syndrome samples using the HaplotypeCaller and GenotypeGVCFS tools for GATK v3.8-0. 
# GATK variant quality score recalibration (VQSR) will be also applied.

use strict;
use warnings;
use File::Basename;

# Please, before running the script first edit the following paths
my $HAPMAP    = "/home/mpinsach/Escriptori/Datasets/hapmap_3.3.hg19.sites.vcf";                           #True variants used for VQSR (type of variants: SNVs)
my $OMNI      = "/home/mpinsach/Escriptori/Datasets/1000G_omni2.5.hg19.sites.vcf";                        #True variants used for VQSR (type of variants: SNVs)
my $HIGHCON   = "/home/mpinsach/Escriptori/Datasets/1000G_phase1.snps.high_confidence.hg19.sites.vcf";    #True variants and false positives used for VQSR (type of variants: SNVs)
my $dbSNP     = "/home/mpinsach/Escriptori/Datasets/dbsnp_138.hg19.vcf";                                  #Variants used for VQSR to stratify output metrics (type of variants: SNVs and Indels)
my $MILLS     = "/home/mpinsach/Escriptori/Datasets/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"; #True variants used for VQSR (type of variants: Indels)
my $reference = "/home/mpinsach/Escriptori/ucsc.hg19/ucsc.hg19.fasta";                                    #Human reference genome (GRCh37/hg19)
my $bed       = "/home/mpinsach/Escriptori/sandbox/RegSeq_numbered_trueregions.bed";                      #List of target regions (n= 1,293 cis-regulatory regions)

my $dir     = $ARGV[0];    #Directory containing aligned BAM files
my $out_vcf = "$dir/vcfs"; #Directory where the final joint VCF will be stored
my $samtools = `which samtools`; 
chomp $samtools;
if (!$samtools) {
    print " ERROR: SAMtools was not found on PATH\n";
    exit;
}
my %GATK = (
	"3.8"=> "/home/mpinsach/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar",
);


if (@ARGV < 2) {
    Help();
}
if (!-e $dir) {
    print " ERROR: Input directory does not exist!\n";
    exit;
}
if (!-e $out_vcf) {
    mkdir $out_vcf;
}

my $varcall_log = "$dir/varcall.log.txt";
open LOG, ">", $varcall_log;

#########################
sub Help {
    print " Usage: $0 <INPUT_DIR> <OUTPUT_DIR>\n";
    exit;
}
#########################

# Globbing inside input directory 

my @bams = glob "$dir/*.bam";
if (!@bams) {
	print " ERROR: no bam files were found\n";
}

my $count = 0;
my $total_samples = scalar (@bams);


open (LOG, ">", $varcall_log) || die "ERROR: unable to open $varcall_log\n";

# Checking bam index
foreach my $bam (@bams) {
	$count++;
	if ( !-e "$bam.bai" ) {
	print " INFO: indexing $bam\n";
		my $cmd = "$samtools index $bam";
		system ($cmd);
	}
	my $name =  basename(( split /\./, $bam )[0]);
	print "Analyzing sample $name ($count of $total_samples)\n";
	my $date = localtime();
	print LOG " INFO: Analyzing sample $name ($count of $total_samples)\n";
	print LOG "$date\n";

# Generating VCF files for each sample
	print "Generating VCF file for sample $name ($count of $total_samples)\n";
		my $cmd = "java -Xmx20g -jar $GATK{3.8} -T HaplotypeCaller -R $reference -L $bed -I $bam --emitRefConfidence GVCF -o $out_vcf/$name.raw.snps.indels.g.vcf";
	system ($cmd);
	
	print LOG "INFO: Generating VCF file for sample $name\n";
	print LOG "$cmd\n";
}

my @vcfs = glob "$out_vcf/*.vcf";
my $str;

foreach my $vcf (@vcfs) {
	$str .= "--variant $vcf ";
}

# Generating a joint variant calling of all samples

print "Running GenotypeGVCF (GATK 3.8-0)\n";
my $cmd = "java -Xmx20g -jar $GATK{3.8} -T GenotypeGVCFs -R $reference $str -o $out_vcf/Regulome.raw.vcf";
system ($cmd);

print LOG " INFO: Running GenotypeGVCF (GATK 3.8-0)\n";
print LOG "$cmd\n";

# Running VariantRecalibrator for SNVs (VQSR)
print "Running VariantRecalibrator for SNVs (VQSR)\n";
$cmd = "java -Xmx20g -jar $GATK{3.8} -T VariantRecalibrator -nt 6 -R $reference -input $out_vcf/Regulome.raw.vcf -L $bed -recalFile $out_vcf/Regulome.snps.recal -tranchesFile $out_vcf/Regulome.snps.tranches -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP -resource:omni,known=false,training=true,truth=false,prior=12.0 $OMNI -resource:1000G,known=false,training=true,truth=false,prior=10.0 $HIGHCON -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbSNP -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff -mode SNP";
system ($cmd);

print LOG " INFO: Running VariantRecalibrator for SNVs (VQSR)\n";
print LOG "$cmd\n";


# Running ApplyRecalibration for SNVs (VQSR)
print "Running ApplyRecalibration for SNVs (VQSR)\n";
$cmd = "java -Xmx20g -jar $GATK{3.8} -T ApplyRecalibration -R $reference -input $out_vcf/Regulome.raw.vcf -recalFile $out_vcf/Regulome.snps.recal -tranchesFile $out_vcf/Regulome.snps.tranches -o $out_vcf/Regulome.recalibrated.snps.vcf -mode SNP";
system($cmd);

print LOG " INFO: Running ApplyRecalibration for SNVs (VQSR)\n";
print LOG "$cmd\n";

# Running VariantRecalibrator for Indels (VQSR)
print "Running VariantRecalibrator for INDELs\n";
$cmd = "java -Xmx20g -jar $GATK{3.8} -T VariantRecalibrator -R $reference -input $out_vcf/Regulome.raw.vcf -recalFile $out_vcf/Regulome.indels.recal -L $bed -tranchesFile $out_vcf/Regulome.indels.tranches --maxGaussians 4 -resource:mills,known=false,training=true,truth=true,prior=12.0 $MILLS -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbSNP -an FS -an SOR -an QD -an ReadPosRankSum -an MQRankSum -an InbreedingCoeff -mode INDEL";
system($cmd);

print LOG " INFO: Running VariantRecalibrator for Indels (VQSR)\n";
print LOG "$cmd\n";

# Running ApplyRecalibration for Indels (VQSR)
print "Running ApplyRecalibration for INDELs (VQSR)\n";
$cmd = "java -Xmx20g -jar $GATK{3.8} -T ApplyRecalibration -R $reference -input $out_vcf/Regulome.recalibrated.snps.vcf -recalFile $out_vcf/Regulome.indels.recal -tranchesFile $out_vcf/Regulome.indels.tranches -o $out_vcf/Regulome.recalibrated.vcf -mode INDEL";
system($cmd);
print LOG " INFO: Running ApplyRecalibration for Indels (VQSR)\n";
print LOG "$cmd\n";
