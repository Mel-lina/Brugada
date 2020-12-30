#!/usr/bin/perl

# Goal: replace 1000 Genomes sample name for Ancestry "_" sample name (ex: AFR_HG01885)
# Script called when running prepare_1KG.sh

use strict;
use warnings;

# Please, before running the script first edit the following paths
my $infile  = "/home/mpinsach/Escriptori/1000Genomes/sequencing/igsr_samples.tsv"; #File containing 1000 Genomes sample names with their respective population descriptions
my $kg_dir  = "/home/mpinsach/Escriptori/1000Genomes/sequencing";                  #Directory where 1000 Genomes VCFs are found
my @chrs    = ("3", "7", "10", "11", "12");                                        #Chromosomes where the targeted 1,293 cis-regulatory regions are found
my %hash    = ();


open (IN, "<", $infile);

while (my $line=<IN>) {
    chomp $infile;

    next if $line =~/^Sample/;
    my @tmp = split (/\t/, $line);

    $hash{$tmp[0]} = $tmp[5];
}
close IN;


foreach my $chr (@chrs) {

    open (IN, "<", "$kg_dir/1000Genomes_chr$chr\_split_snvs_target_norm.vcf") || die " ERROR: $kg_dir/1000Genomes_chr$chr\_split_snvs_target_norm.vcf does not exist\n";
    open (OUT_KG, ">", "$kg_dir/1000Genomes_chr$chr\_split_snvs_target_norm_pop.vcf");
    while (my $line=<IN>) {
        chomp $line;

    if ($line=~/^##/) {
        print OUT_KG "$line\n";
        next;
    }
    if ($line =~/^#CHROM/) {
        my @tmp = split (/\t/, $line);
        my @tmpHeader = ();
        for (my $i=9; $i<@tmp; $i++) {
            my $newName = $hash{$tmp[$i]} . "_" . $tmp[$i];
            push @tmpHeader, $newName;
        }
        my $header = join ("\t", @tmpHeader);
        my $befSampleNames = join ("\t", @tmp[0..8]);
        print OUT_KG "$befSampleNames\t$header\n";
        next;
    }

    print OUT_KG "$line\n";
}
}
