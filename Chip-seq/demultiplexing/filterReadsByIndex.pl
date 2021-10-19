#!/usr/bin/perl -w

use strict;
use warnings;


use Getopt::Long;
use Pod::Usage;

# perl /home/hiseq/scripts/demultiplexing/filterReadsByIndex.pl --indexToSkip index.txt --1 Ara-B1A1-2_S4_L001_R1_001.fastq.gz --2 Ara-B1A1-2_S4_L001_R2_001.fastq.gz

my %options;

GetOptions(\%options, 'indexToSkip=s', '1=s', '2=s', 'help|h') or die("Error in command line arguments\n");


if($options{'help'} || $options{'h'}){
	&pod2usage({EXIT => 2, VERBOSE => 2});
	#exit 1;
}


if(!$options{'indexToSkip'}){
	print STDERR "Error: Please provide the index file\n";
	&pod2usage({EXIT => 2, VERBOSE => 0});
}


if(!$options{'1'}){
	print STDERR "Error: Please provide the Read_1 file\n";
	&pod2usage({EXIT => 2, VERBOSE => 0});
}


if(!$options{'2'}){
	print STDERR "Error: Please provide the Read_2 file\n";
	&pod2usage({EXIT => 2, VERBOSE => 0});
}


my $isGz = 0;
if($options{1}=~m/gz$/){
	$isGz=1;
}


open(my $idxFh, $options{'indexToSkip'}) or die "Cannot open file $options{'indexToSkip'}: $!";

my $pattern = join("|", <$idxFh>);

#decide new file names and create them
$options{1}=~m/(.*)\.fastq.*/;
my $newR1 = $1."_indexFiltered.fastq.gz";
$options{2}=~m/(.*)\.fastq.*/;
my $newR2 = $1."_indexFiltered.fastq.gz";

open(my $fhR1, "|gzip >$newR1") or die "Cannot create file $newR1: $!";
open(my $fhR2, "|gzip >$newR2") or die "Cannot create file $newR2: $!";


#Read fastq data
open(my $fh1, $isGz ? "gzip -dc $options{1} |" : $options{1}) or die "Cannot open file $options{1}: $!";
open(my $fh2, $isGz ? "gzip -dc $options{2} |" : $options{2}) or die "Cannot open file $options{2}: $!";

my $count = 0;
my ($p1, $p2);

while(1){
	if(eof($fh1)){
		last;
	}
	
		#read line1: headers
	$p1 = <$fh1>;
	$p2 = <$fh2>;
	
	if($p1=~m/^@\w+:\d+:[\w-]+(:\d+){4}\s\d+:(Y|N):0:/ && $p2=~m/^@\w+:\d+:[\w-]+(:\d+){4}\s\d+:(Y|N):0/){
		#@D00691:42:HCVYCBCXX:1:1107:1312:2202 1:N:0:TTAGGC
		#read sequence
		if($p1 !~ /\s\d+:(Y|N):0:($pattern)/){
			$p1 .= <$fh1>;
			$p2 .= <$fh2>;
			
			#read + line
			$p1 .= <$fh1>;
			$p2 .= <$fh2>;
			
			#read qual line
			$p1 .= <$fh1>;
			$p2 .= <$fh2>;
			
			print $fhR1 $p1;
			print $fhR2 $p2;
			
			$count++;
		}
	}
}


close($fhR1);
close($fhR2);
close($fh1);
close($fh2);
close($idxFh);

