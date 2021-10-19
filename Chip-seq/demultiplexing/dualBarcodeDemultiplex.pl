#!/usr/bin/perl -w

use strict;
use warnings;


use Getopt::Long;
use Pod::Usage;
use File::Basename;
use List::MoreUtils qw(uniq);


# perl /home/lakhanp/scripts/demultiplexing/dualBarcodeDemultiplex.pl --barcodes FHS_mix28DualBarcodes.txt -1 rawData/FHSmix28_S1_L001_R1_001.fastq.gz -2 rawData/FHSmix28_S1_L001_R2_001.fastq.gz

my %options;
my $isPaired = 1;

GetOptions(\%options, 'barcodes=s', '1=s', '2=s', 'help|h') or die("Error in command line arguments\n");


if($options{'help'} || $options{'h'}){
	&pod2usage({EXIT => 2, VERBOSE => 2});
	#exit 1;
}


if(!$options{'barcodes'}){
	print STDERR "Error: Please provide the barcode file\n";
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


#read barcode data:
my %barcodes = ();
my %barLength = ();
my %files = ();

my @r1Bars = ();
my @r2Bars = ();

open(my $barFh, $options{'barcodes'}) or die "Cannot open file $options{'barcodes'}: $!";


while(<$barFh>){
	if(/^(\w\S+)\s+(\w+)\s+(\w+)\s*$/){
		my $name = $1;
		my $bar1 = $2;
		my $bar2 = $3;
		
		$bar1 = uc $bar1;
		$bar2 = uc $bar2;
		
		push(@r1Bars, $bar1);
		push(@r2Bars, $bar2);
		
		my $bar = $bar1.'_'.$bar2;
				
		#check if the same file name is used for different barcodes
		if(exists $files{$name}){
			print "Error: Two different barcodes are using same sample name.
			$files{$name} : $name
			$bar : $name\n";
			die;
		}
		
		$files{$name} = $bar;
		
		#Check if the barcode pair is repeated
		if(exists $barcodes{$bar}){
			print "Error: Two different samples have same pair of barcodes
			Sample $barcodes{$bar}->{id} : $barcodes{$bar}->{bar1}, $barcodes{$bar}->{bar2}
			Sample $name : $bar1, $bar2\n";
			die;
		}
		
		
		$barcodes{$bar}->{'id'} = $name;			#Sample name prefix
		$barcodes{$bar}->{'bar1'} = $bar1;			#Sample name prefix
		$barcodes{$bar}->{'bar2'} = $bar2;			#Sample name prefix
		
		open($barcodes{$bar}->{'fhR1'}, "|gzip >$barcodes{$bar}->{id}_R1.fastq.gz") or die "Cannot create file $barcodes{$bar}->{id}_R1.fastq.gz: $!";			#R1 file
		open($barcodes{$bar}->{'fhR2'}, "|gzip >$barcodes{$bar}->{id}_R2.fastq.gz") or die "Cannot create file $barcodes{$bar}->{id}_R2.fastq.gz: $!";			#R2 file		
		
		#this barcode length will be used to trim the sequence and qual line
		$barLength{$bar1} = length($bar1);
		$barLength{$bar2} = length($bar2);
		
		$barcodes{$bar}->{'count'} = 0;			#counter for each barcode
	}
	elsif(/^\s*$/){
		next;
	}
	else{
		print "Wrong format is barcode file at line: $_";
		die;
	}
}


# $options{1}=~m/(.*)\.fastq.*/;
my $unknownR1 = basename($options{1});
$unknownR1 = 'unknown_'.$unknownR1;

my $unknownR2 = basename($options{2});
$unknownR2 = 'unknown_'.$unknownR2;

$barcodes{'unknown'}->{'id'} = 'unknown';
open($barcodes{'unknown'}->{'fhR1'}, "|gzip >$unknownR1") or die "Cannot create file $unknownR1: $!";
open($barcodes{'unknown'}->{'fhR2'}, "|gzip >$unknownR2") or die "Cannot create file $unknownR2: $!";


$barLength{'unknown'} = 0;				#barcode length for unknown sequences
$barcodes{'unknown'}->{'count'} = 0;				#counter for unknown sequences

my $pattern1 = join("|", uniq(@r1Bars));
my $pattern2 = join("|", uniq(@r2Bars));

print "Barcode1: $pattern1\nBarcode2: $pattern2\n";

#Read fastq data
open(my $fh1, $isGz ? "gzip -dc $options{1} |" : $options{1}) or die "Cannot open file $options{1}: $!";
open(my $fh2, $isGz ? "gzip -dc $options{2} |" : $options{2}) or die "Cannot open file $options{2}: $!";



my ($p1, $p2, $bar1, $bar2, $outBar);

while(1){
	if(eof($fh1)){
		last;
	}
	
	#read line1: headers
	$p1 = <$fh1>;
	$p2 = <$fh2>;
	
	if($p1=~m/^@\w+:\d+:[\w-]+(:\d+){4}\s\d+:(Y|N)/ && $p2=~m/^@\w+:\d+:[\w-]+(:\d+){4}\s\d+:(Y|N)/){
		#@SIM:1:FCX:1:15:6329:1045 1:N:0:2
		
		#read line 2: sequence
		my $sq1 = <$fh1>;
		my $sq2 = <$fh2>;
		
		if($sq1=~m/^($pattern1)/){
			$bar1 = $1;
			
			#if bar1 is found at 5' end of R1 AND bar2 is found at 5' end of R2 then only go ahead.
			if($sq2=~m/^($pattern2)/){
				$bar2 = $1;
				
				$outBar = $bar1.'_'.$bar2;
				
				#all combinations of bar1_bar2 will not be used in multiplexing. So ignore the bar1_bar2 combination for which there is no sample ID. Process it as unknown
				if(!exists $barcodes{$outBar}){
					#print $barcodes{$outBar}->{'id'},"\t$bar1\t$bar2\t$outBar\n";
					#print $p1,$sq1,$p2,$sq2;
					$bar1 = 'unknown';
					$bar2 = 'unknown';
					$outBar = 'unknown';
					$p1 .= $sq1;
					$p2 .= $sq2;
				}
				else{
				#trim from 5_prime end to remove barcode
					$p1 .= substr($sq1, $barLength{$bar1});
					$p2 .= substr($sq2, $barLength{$bar2});
				}
			}
			else{
				$bar1 = 'unknown';
				$bar2 = 'unknown';
				$p1 .= $sq1;
				$p2 .= $sq2;
				$outBar = 'unknown';
			}
		}
		else{
			$bar1 = 'unknown';
			$bar2 = 'unknown';
			$p1 .= $sq1;
			$p2 .= $sq2;
			$outBar = 'unknown';
		}
		
		#read line 3: +
		$p1 .= <$fh1>;
		$p2 .= <$fh2>;
		
		#read line 4: qual
		#trim the qual of 5_prime barcode
		$p1 .= substr(<$fh1>, $barLength{$bar1});
		$p2 .= substr(<$fh2>, $barLength{$bar2});
		
		print {$barcodes{$outBar}->{'fhR1'}} $p1;
		print {$barcodes{$outBar}->{'fhR2'}} $p2;
		$barcodes{$outBar}->{'count'}++;
	}
}




foreach(sort{$barcodes{$a}->{'id'} cmp $barcodes{$b}->{'id'}}keys %barcodes){
	close($barcodes{$_}->{'fhR1'});
	close($barcodes{$_}->{'fhR2'});
	
	print "$barcodes{$_}->{id}\t$barcodes{$_}->{count}\n";
}


close($barFh);







__END__


=head1 NAME


=head1 SYNOPSIS

perl dualBarcodeDemultiplex.pl --barcodes <Barcode file> -1 <Mate1> -2 <Mate2>

Help Options:

	--help	Show this scripts help information.

=head1 DESCRIPTION

This script demultiplex the dual barcoded paired end reads. Please provide the sample Id and barcode pair information as tab separated file.
E.g.: 	ID1	ATGGC	GGCAT
		
		ID2	GGCAT	GTTAT
		
		.	.		.

=head1 OPTIONS

=over 30

=item B<--barcodes>

[STR] Barcode file

=item B<--1>

[STR] Forward read file

=item B<--2>

[STR] Reverse read file

=item B<--help>

Show this scripts help information.

=back


=cut




