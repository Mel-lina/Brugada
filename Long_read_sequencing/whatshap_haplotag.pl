#! /usr/bin/perl

# Goal: tag phased haplotypes for visualization
# Haplotagging is performed using WhatsHap.

use strict;
use warnings;
use diagnostics;
use File::Basename;

# Please, before running the script first edit the following paths
my $in_dir = $ARGV[0];	                                                #Directory containing filtered BAMs obtained with the script MinION_alignment_and_filtering.pl
my $out_dir = $ARGV[1];                                                 #Directory where WhatsHap haplotagged VCF files will be stored
my $REF_GENOME = "/home/mpinsach/Escriptori/ucsc.hg19/ucsc.hg19.fasta"; #Reference genome (GRCh37/hg19)


if (@ARGV < 3) {
    Help();
}
if (!-e $in_dir) {
    print " ERROR: Input directory does not exist!\n";
    exit;
}
if (!-e $out_dir) {
    mkdir $out_dir;
}

my $logfile = "$output_dir/process_log.txt";
open LOG, ">", $logfile;

sub Help {
    print " Usage: $0 <INPUT_DIR> <OUTPUT_DIR> <REF_GENOME>\n";
    exit;
}

my $whatshap = `which whatshap`; 
chomp $whatshap;
if (!$whatshap) {
    print " ERROR: whatshap was not found on PATH\n";
    exit;
}

my $samtools = `which samtools`; 
chomp $samtools;
if (!$samtools) {
    print " ERROR: SAMtools was not found on PATH\n";
    exit;
}

# Globing inside input directory 

my @bams = glob "$in_dir/*.bam";
if (!@bams) {
	print " ERROR: no bam files were found\n";
}

my $count = 0;
my $total_samples = scalar (@bams);


# Checking BAM index
foreach my $bam (@bams) {
	$count++;
	if ( !-e "$bam.bai" ) {
	print " INFO: indexing $bam\n";
		my $cmd = "$samtools index $bam";
		system ($cmd);
	}
	my $bam_name = basename($bam);
	my @tmp = split("_", $bam_name);
	my $name = $tmp[0];
	my $barcode     = $tmp[1];
	
# Running WhatsHap haplotag to get haplotagged VCFs
	print "INFO: Running WhatsHap haplotag for sample $name\n";
	my $cmd1 = "$whatshap haplotag -o $out_dir/$name\_$barcode\_haplotagged_no_clipping_length8kb.bam --reference $REF_GENOME  /home/mpinsach/Escriptori/MinION/whatshap/phase/$name\_$barcode\_phased.vcf.gz $bam --ignore-read-groups";
	system ($cmd1);
	print LOG "INFO: Running WhatsHap haplotag for sample barcode_$i\n";
    print LOG " INFO: $cmd1\n";
}

