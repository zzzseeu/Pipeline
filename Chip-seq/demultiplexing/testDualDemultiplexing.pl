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

GetOptions(\%options, 'barcode=s', '1=s', '2=s', 'out=s', 'help|h') or die("Error in command line arguments\n");


if($options{'help'} || $options{'h'}){
	&pod2usage({EXIT => 2, VERBOSE => 2});
	#exit 1;
}


if(!$options{'barcode'}){
	print STDERR "Error: Please provide the barcode\n";
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

if(!$options{'out'}){
	print STDERR "Error: Please provide the output prefix\n";
	&pod2usage({EXIT => 2, VERBOSE => 0});
}


my $isGz = 0;
if($options{1}=~m/gz$/){
	$isGz=1;
}



my $pattern1 = $options{'barcode'};


#Read fastq data
open(my $fh1, $isGz ? "gzip -dc $options{1} |" : $options{1}) or die "Cannot open file $options{1}: $!";
open(my $fh2, $isGz ? "gzip -dc $options{2} |" : $options{2}) or die "Cannot open file $options{2}: $!";

my $r1 = $options{'out'}.'_R1.fastq.gz';
my $r2 = $options{'out'}.'_R2.fastq.gz';

open(my $out1, "|gzip >$r1") or die "Cannot create file $r1: $!";
open(my $out2, "|gzip >$r2") or die "Cannot create file $r2: $!";

my ($p1, $p2);

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
			$p1 .= $sq1;
			$p2 .= $sq2;
			
			#read line 3: +
			$p1 .= <$fh1>;
			$p2 .= <$fh2>;
			
			#read line 4: qual
			#trim the qual of 5_prime barcode
			$p1 .= <$fh1>;
			$p2 .= <$fh2>;
			
			print $out1 $p1;
			print $out2 $p2;
		}
		else{
			$p1 .= $sq1;
			$p2 .= $sq2;
			
			#read line 3: +
			$p1 .= <$fh1>;
			$p2 .= <$fh2>;
			
			#read line 4: qual
			#trim the qual of 5_prime barcode
			$p1 .= <$fh1>;
			$p2 .= <$fh2>;
		}
	
	}
}




close($fh1);
close($fh2);
close($out1);
close($out2);








__END__


=head1 NAME


=head1 SYNOPSIS

perl dualBarcodeDemultiplex.pl --barcodes <Barcode file> -1 <Mate1> -2 <Mate2> --out <Out file prefix>

Help Options:

	--help	Show this scripts help information.

=head1 DESCRIPTION

This script is specifically written to test the dual barcode demultiplexing script. It takes on barcode sequence and extracts all the reads that begins with the barcode sequence from R1 file. It writes the respective R2 reads in outPrefix_R2 file.
		.	.		.

=head1 OPTIONS

=over 30

=item B<--barcodes>

[STR] Barcode file

=item B<--1>

[STR] Forward read file

=item B<--2>

[STR] Reverse read file

=item B<--out>

[STR] Output file prefix

=item B<--help>

Show this scripts help information.

=back


=cut




