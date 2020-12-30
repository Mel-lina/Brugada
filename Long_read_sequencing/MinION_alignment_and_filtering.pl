#!/usr/bin/perl

# Goal: alignment and filtering of long-read sequencing data (generated on the MinION)
# Alignment is performed using minimap2, which requires a minimap2 index before runing the algorithm
# Filtering consists in the removal of: (1) soft and hard-clipped reads, (2) supplementary aligments, and (3) reads smaller than 8kb

use strict;
use warnings;
use diagnostics;
use File::Basename;

# Please, before running the script first edit the following paths
my $input_dir     = $ARGV[0];   #Directory containing MinION Fastq files stored in different folters, depending on the sequencing barcode used
my $output_dir    = $ARGV[1];   #Directory where aligned and filter BAM files will be sotred
my $minimap2_index = $ARGV[2];  #Before runing the script, you must create a minimap2 index

if (@ARGV < 3) {
    Help();
}
if (!-e $input_dir) {
    print " ERROR: Input directory does not exist!\n";
    exit;
}
if (!-e $minimap2_index) {
    print " ERROR: Minimap2 index was not found\n";
}
if (!-e $output_dir) {
    mkdir $output_dir;
}

my $logfile = "$output_dir/process_log.txt";
open LOG, ">", $logfile;

sub Help {
    print " Usage: $0 <INPUT_DIR> <OUTPUT_DIR> <MINIMAP_INDEX_FILE>\n";
    exit;
}

my $minimap2 = `which minimap2`; 
chomp $minimap2;
if (!$minimap2) {
    print " ERROR: minimap2 was not found on PATH\n";
    exit;
}

my $samtools = `which samtools`; 
chomp $samtools;
if (!$samtools) {
    print " ERROR: SAMtools was not found on PATH\n";
    exit;
}

# Scheme required: input_dir/barcode1/*.fastq.gz

my @dirs = glob ("$input_dir/*");
@dirs = grep ( -d $_, @dirs );

if (!@dirs) {
    print " ERROR: No sample directories were found at $input_dir\n";
    exit;
}

my $i = 0;
foreach my $dir (@dirs) {
    $i++;
    my @fastq = glob ("$dir/*.fastq");
    if (!@fastq) {
        print " Skipping directory $dir\n";
        next;
    }
    print " $dir\n";

    my $fastq1 = basename($fastq[0]);
    my @tmp = split("_", $fastq1);
    my $barcode = $tmp[2];

    my $fastq_name = basename($fastq);
	my @tmp = split("_", $fastq_name);
	my $name = $tmp[0];
	my $barcode = $tmp[1];

    # Concatenate all fastq files for the same sample into a single gzipped fastq file
    print "INFO: Concatenating raw fastq into a single gzipped fastq file for sample $name\n";
    my $cmd = "cat $dir/*.fastq | gzip - > $dir/$name\_$barcode.fastq.gz";
    system $cmd;
    print LOG " INFO: Concatenating raw fastq into a single gzipped fastq file for sample $name\n";
    print LOG " INFO: $cmd\n";

    # Align long reads with minimap2
    $cmd = "minimap2 -x map-ont -t6 -a $minimap2_index $dir/$name\_$barcode.fastq.gz | samtools sort -O BAM - -o $output_dir/$name\_$barcode.bam";
    system $cmd;
    print LOG " INFO: Minimap2 alignment for sample $output_dir/$name\_$barcode.bam\n";
    print LOG " INFO: $cmd\n";

    # Index aligned BAM files
    $cmd = "samtools index $output_dir/$name\_$barcode.bam";
    system $cmd;
    print LOG " INFO: Indexint $output_dir/$name\_$barcode.bam\n";
    print LOG " INFO: $cmd\n";

    # Remove soft and hard-clipped reads
    print " INFO: Removing soft/hard-clipped reads\n";
    print LOG " INFO: Removing soft/hard-clipped reads\n";
    $cmd = "samtools view $output_dir/$name\_$barcode.bam | awk \'{if (\$6 !~/S|H/){print \$0} }\' > $output_dir/$name\_$barcode.body.sam";
    system $cmd;

    print " INFO: Extracting header for $name bam\n";
    print LOG " INFO: Extracting header for $name bam\n";
    $cmd = "samtools view -H $output_dir/$name\_$barcode.bam > $output_dir/$name\_$barcode.header.sam\n";
    system $cmd;
    print LOG " INFO: $cmd\n";

    $cmd = "cat $output_dir/$name\_$barcode.header.sam $output_dir/$name\_$barcode.body.sam > $output_dir/$name\_$barcode.sam";
    system $cmd;
    print LOG " INFO: $cmd\n"; 
    unlink ("$output_dir/$name\_$barcode.header.sam", "$output_dir/$name\_$barcode.body.sam");

    # Remove reads <8kb and supplementary aligments
    print " INFO: Removing short reads (<8 kb) and supplementary/secondary aligments\n";
    print LOG "INFO: Removing short reads (<8 kb) and supplementary/secondary aligments\n";
    $cmd = "samtools view -bhS $output_dir/$name\_$barcode.sam| samtools view -bh -F 2304 | samtools view -h | awk 'length(\$10) > 8000 || \$1 ~ /^@/' | samtools view -bS - > $output_dir/$name\_$barcode\_no_cliping_length8kb.bam";
    system $cmd;
    print LOG " INFO: $cmd\n";

    $cmd = "samtools index $output_dir/$name\_$barcode\_no_cliping_length8kb.bam";
    system $cmd;
    print LOG " INFO: $cmd\n";

}
