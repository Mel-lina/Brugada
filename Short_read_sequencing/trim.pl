#!/usr/bin/perl

# Goal: Trimming of short reads (called in the Short_read_alignmenmt.pl script)

use strict;
use warnings;

my $infile=$ARGV[0]; #Directory containing FASTQ files 
my $outfile="$infile";
$outfile=~s/fastq/trimmed.fastq/;
$outfile=~s/.*\///;
my $logfile=$outfile;
$logfile=~s/fastq/log/;
if (!-e $infile){
	print "##ERROR\tINFILE $infile not found\n";
	exit;
}
if ($infile =~ /\.gz$/) {
	open(IN, "gunzip -c $infile |") || die "can't open pipe to $infile";
	$outfile=~s/.gz//;
	$logfile=~s/.gz//;
}
else {
	open(IN, $infile) || die "can't open $infile";
}

open OUT, ">$outfile";
my $count=0;
my $header='';
my $seq='';
my $qual='';
my $joint='';
my $readbase=0;
my $trimbase=0;
my $finalbase=0;
my $reads=0;
my $readlength=0;
while(my $inl=<IN>){
	if ($count==0)
	{
		$header=$inl;	
	}
	elsif ($count==1)
	{
		chomp $inl;
		$seq=$inl;
	}
	elsif ($count==2)
	{
		$joint=$inl;
	}	
	else{
		$reads++;
		$count=0;
		chomp $inl;
		$qual=$inl;
		my $length=length $qual;
		$readlength=$length;
		$readbase+=$length;
		my @tmp=split //, $qual;
		my $i=0;
		for ($i=$length-1; $i>=0; $i--){
			last if(ord $tmp[$i]>43);
		}
		if ($i<39)
		{
			$i=39;
		}
		$i++;
		$trimbase+=$length-$i;
		$finalbase+=$i;
		my $finseq=substr $seq,0, $i;
		my $finqual=substr $qual,0,$i;
		print OUT "$header$finseq\n$joint$finqual\n";
		next;
	}
	$count++;
}
close OUT;
close IN;
open  LOG, ">$logfile" || die "cannot open $logfile file\n";
my $tmp=0;
print LOG "readlength;$readlength\n";
print LOG "readbase;$readbase\n";
if ($readbase==0){
	$tmp=sprintf "%6.2f" ,0;
}
else
{
	$tmp=sprintf "%6.2f" ,100*$trimbase/$readbase;
}
print LOG "trimbase;$trimbase;$tmp;\n";
if ($readbase==0){
	$tmp=sprintf "%6.2f" ,0;
}
else
{
	$tmp=sprintf "%6.2f" ,100*$finalbase/$readbase;
}
print LOG "finalbase;$finalbase;$tmp;\n";
print LOG "reads;$reads\n";
close LOG;

