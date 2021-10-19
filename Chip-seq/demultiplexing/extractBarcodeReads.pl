#!/usr/bin/perl -w

use strict;
use warnings;

if($#ARGV == -1){
	print STDERR "Please provide input FASTQ file\n";
	die
}

my %barcodes = ();
while(<DATA>){
	chomp;
	my @tmp = split(/\t/, $_);
	$barcodes{$tmp[1]}->{'file'} = $tmp[0].'_'.$tmp[1].'.fastq.gz';
	open($barcodes{$tmp[1]}->{'fh'}, '|gzip >'.$barcodes{$tmp[1]}->{'file'}) or die "Cannot open file ",$barcodes{$tmp[1]}->{'file'},": $!\n";
}

my $pattern = join('|', keys %barcodes);

open(my $fh, $ARGV[0]) or die "Cannot open file $ARGV[0]: $!\n";
my $p1 = '';

while(1){
	if(eof($fh)){
		last;
	}
	
	$p1 = <$fh>;
	if($p1=~m/^\@HWI/){
		if($p1=~m/:\d+#($pattern)\//){
			$p1 .= <$fh>;
			$p1 .= <$fh>;
			$p1 .= <$fh>;
			
			print {$barcodes{$1}->{'fh'}} $p1;
			$barcodes{$1}->{'count'}++;
			#print $1,"\n";
		}
	}
}


foreach(keys %barcodes){
	print $_,"\t",$barcodes{$_}->{'count'},"\n";
}


close($fh);









__DATA__
bar1	GGCTAC
bar2	CTTGTA
