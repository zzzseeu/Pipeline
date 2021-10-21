#!/usr/bin/perl -w

use strict;
use warnings;

if(@ARGV != 1){
	die "ERROR: Please provide fastq file name\n";
}

open(my $fh, $ARGV[0]=~m/gz$/ ? "gzip -dc $ARGV[0] |" : $ARGV[0]) or die "Cannot open file $ARGV[0]: $!";

my %enrichment = ();
my $i = 3;
while(<$fh>){
	
	if($i%4==0){
		/^(\w{8})/;
		$enrichment{$1}++;
	}
	
	$i++;
}

close($fh);


foreach(sort{$enrichment{$b} <=> $enrichment{$a}}keys %enrichment){
	print "$_\t$enrichment{$_}\n";
}


