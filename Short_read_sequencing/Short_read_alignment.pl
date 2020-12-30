#!/usr/bin/env perl

# Goal: Preprocess FASTQ files from all 86 Brugada syndrome samples and align them to the human reference genome GRCh37/hg19.

use strict;
use warnings;
use diagnostics;
use Getopt::Long;
use File::Basename;
use Cwd;
use Sort::Key::Natural qw(natsort);



# Please, before running the script first edit the following paths
my $nextera_adaptors = "~/skewer.0.1.123/nextera_adaptors.fa";                     #Adaptors used during library preparation (Nextera Rapid Capture Custom Enrichment Kit)
my $trimBases        = "trim.pl";                                                  #In-house perl script used to quality-prune read ends
my $reference        = "/home/mpinsach/Escriptori/ucsc.hg19/ucsc.hg19.fasta";      #Reference genome (GRCh37/hg19)
my $skewer           = "~/skewer.0.1.123/skewer-0.1.123-linux-x86_64";             #Algorithm used to remove adaptor sequences added during library preparation
my $picard           = "~/Picard/picard.jar";                                      #Algorithm used to remove sequencing duplicates
my $samtools         = "";                                                         #Algorithm used to remove reads with ambiguous multiple secondary alignments
my $gatk             = "/home/ubuntu/GenomeAnalysisTK-3.8-0/GenomeAnalysisTK.jar"; #Agorithm used to left-align indels
my $bwa              = "";                                                         #Algorithm used to align sequencing reads to the human reference genome GRCh37/hg19 


my $fastq_dir; #Directory containing raw FASTQ files 
my $out_dir;   #Directory where aligned BAM files will be stored

if (@ARGV < 2) {
    Help();
}

 # Check input and output directories/files
 if (!$fastq_dir) {
    print " ERROR: no input FASTQ directory was introduced\n";
    Help();
 }
 if (!-e $fastq_dir) {
    print " ERROR: FASTQ directory introduced does not exist!\n";
    Help()();
 }
 if (!$reference) {
    print "ERROR: no reference genome was introduced\n";
    Help()();
 }
 if (!$out_dir ) {
	print "ERROR: no output directory was introduced\n";
    Help()();
 }
 if (!-e $out_dir ) {
	mkdir $out_dir;
 }

 ####################
 sub Help {
    print " Usage: $0 <INPUT_DIR> <OUTPUT_DIR>\n";
    exit;
}
####################

my $align_log = "$out_dir/align.log.txt";
open LOG, ">", $align_log;

######## Bulk processing

 # Get all FASTQ files from fastq_dir
 my @FASTQ = glob ("$fastq_dir/*fastq.gz");
 if (!@FASTQ) {
     print " ERROR: No FASTQ files found at $fastq_dir\n";
     exit;
 }
 
 # Step 1. Adapter trimming with skewer
 my %seenFastq = ();
 foreach my $fastq (@FASTQ) {
     my $name = basename($fastq);
     my @tmp = split ('_', $name);
     $name = $tmp[0];
     
    $seenFastq{$name}++;
    next if $seenFastq{$name} > 1;

     my $fq1; 
     my $fq2;
     if ($fastq =~/R1/) {
         $fq1 = $fastq;
         $fq2 = $fq1;
         $fq2 =~s/R1/R2/;
     }  
     else {
         $fq2 = $fastq;
         $fq1 = $fq2;
         $fq1 =~s/R2/R1/;   
     }

    my $cmd = " -x $nextera_adaptors -m pe $fq1 $fq2 -o $out_dir/$name -t 6 1 > $out_dir/$name.skewer.log";
    system $cmd;
    print LOG " INFO: Adapter trimming of sample $name with Skewer\n";
    print LOG "$cmd\n";

    my $fq1_skewer = "$out_dir/$name-trimmed-pair1.fastq";
    my $fq2_skewer = "$out_dir/$name-trimmed-pair2.fastq";

    # Step 2. Trimming of low quality bases
    $cmd = "perl $trimBases $fq1_skewer"; # For fq1
    system $cmd;
    print LOG " INFO: Trimming low quality bases for $fq1 with trim.pl\n";
    print LOG "$cmd\n";

    $cmd = "perl $trimBases $fq2_skewer"; # For fq2
    system $cmd;
    print LOG " INFO: Trimming low quality bases for $fq2 with trim.pl\n";
    print LOG "$cmd\n";

    unlink $fq1_skewer, $fq2_skewer;

    my $old_fq1_trimmed = "$out_dir/$name-trimmed-pair1.trimmed.fastq";
    my $fq1_trimmed = "$out_dir/$name.trimmed.R1.fastq";
    rename $old_fq1_trimmed, $fq1_trimmed;

    my $old_fq2_trimmed = "$out_dir/$name-trimmed-pair2.trimmed.fastq";
    my $fq2_trimmed = "$out_dir/$name.trimmed.R2.fastq";
    rename $old_fq2_trimmed, $fq2_trimmed;

    # Step 3. Mapping with BWA mem
	$cmd = "$bwa mem $reference -M -t 4 $fq1_trimmed $fq2_trimmed | $samtools view  -Shu -F 256 - | $samtools sort -T $name -o $out_dir/$name.bam";
	system $cmd;
    print LOG " INFO: Read alignment for $name using BWA-MEM\n";
    print LOG "$cmd\n";

    $cmd = "$samtools index $out_dir/$name.bam";
    system $cmd;
  
    # Step 4. Removal of PCR duplicates
    $cmd = "java -jar $picard MarkDuplicates I=$out_dir/$name.bam O=$out_dir/$name.nodup.bam M=$out_dir/$name.picard.metrics.log REMOVE_DUPLICATES=true ASSUME_SORTED=true";
    system $cmd;
    print LOG " INFO: Removing sequencing duplicates for $name using Picard\n";
    print LOG "$cmd\n";

    $cmd = "$samtools index $out_dir/$name.nodup.bam";
    system $cmd;

    # Step 5. Left align INDELS
    $cmd =  "java -Xmx20g -jar $gatk -T LeftAlignIndels -R $reference -I $out_dir/$name.nodup.bam -o $out_dir/$name.nodup.leftAlindel.bam";
    system $cmd;
    print LOG " INFO: Indel left-alignment for $name using GATK v.3.8-0\n";
    print LOG "$cmd\n";

    $cmd = "$samtools index $out_dir/$name.nodup.leftAlindel.bam";
    system $cmd;

 }