#!/usr/bin/perl

# Goal: replace Brugada syndrome sample name for Ancestry "_" sample name (ex: BrS_Sample1)
# Script called when running ancestry_BrS.sh

use strict;
use warnings;

# Please, before running the script first edit the following paths
my $infile = "/home/mpinsach/Escriptori/Ancestry/BrS/BrS_population.txt"; #File containing Brugada syndrome sample names with their respective population description (i.e. BrS)
my $vcf    = "/home/mpinsach/Escriptori/Ancestry/BrS/BrS_86_split_SNVs_noasterisk_norm.vcf"; #Brugada syndrome VCF file generated using the script ancestry_BrS.sh

my %hash = ();

open (IN, "<", $infile);

while (my $line=<IN>) {
    chomp $infile;

    next if $line =~/^Sample/;
    my @tmp = split (/\t/, $line);
    $hash{$tmp[0]} = $tmp[1];
}
close IN;

open (VCF, "<", $vcf);
while (my $line=<VCF>) {
    chomp $line;

    if ($line=~/^##/) {
        print "$line\n";
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
        print "$befSampleNames\t$header\n";
        next;
    }

    print "$line\n";

}
