#!/usr/bin/env perl 

use strict;
use warnings;


use Getopt::Long;
use Pod::Usage;
use File::Basename;



my %options;
$options{'suffix'} = '';
my $isPaired = 1;
$options{'barTrim'} = 7;

GetOptions(\%options, 'barcodes=s', 'barTrim=i', '1=s@', '2:s@', 'suffix=s', 'help|h') or die("Error in command line arguments\n");


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
	print STDERR "Warning: Read 2 file not provided. Assuming single end data\n";
	$isPaired = 0;
}

#############################################################################
## validation
my %barcodes = ();
my %files = ();
my $barTrim = $options{'barTrim'};

my @r1Files = split(/,/, join(',', @{$options{1}}));
my @r2Files = $isPaired ? split(/,/, join(',', @{$options{2}})) : ();

print STDERR "R1 files: @r1Files\nR2files: @r2Files\n";

my $isGz = 0;
foreach(@r1Files, @r2Files){
	# print $_,"\n";
	if($_!~m/gz$/){
		die "ERROR: Only .gz files accepted. File $_ is not a gzip compressed file.\n";
	}
	
	if(! -e $_){
		die "ERROR: missing file $_\n";
	}
}

if($isPaired && $#r1Files != $#r2Files){
	die "Error: Number of R1 and R2 files not same\n";
}

## undetermined reads files
$barcodes{'unknown'}->{'name'} = 'unknown'.'_'.$options{'suffix'};
$barcodes{'unknown'}->{'R1'} = $options{'suffix'} eq '' ? 'unknown_R1.fastq.gz' : 'unknown_'.$options{'suffix'}.'_R1.fastq.gz';
$barcodes{'unknown'}->{'R2'} = $options{'suffix'} eq '' ? 'unknown_R2.fastq.gz' : 'unknown_'.$options{'suffix'}.'_R2.fastq.gz';


open($barcodes{'unknown'}->{'fhR1'}, "|gzip > $barcodes{'unknown'}->{'R1'}") or die "Cannot create file $barcodes{'unknown'}->{'R1'}: $!";
$barcodes{'unknown'}->{'barLen'} = 0;				#barcode length for unknown sequences
$barcodes{'unknown'}->{'count'} = 0;				#counter for unknown sequences

#Read fastq data
my $r1OpenCmd = join(' ', 'gzip -dc', @r1Files, '|');
my $r2OpenCmd = join(' ', 'gzip -dc', @r2Files, '|');

# print "$r1OpenCmd\n$r2OpenCmd\n";

open(my $fh1, $r1OpenCmd) or die "Cannot open file/s @r1Files: $!";
my $fh2 = undef;

if($isPaired){
	open($barcodes{'unknown'}->{'fhR2'}, "|gzip > $barcodes{'unknown'}->{'R2'}") or die "Cannot create file $barcodes{'unknown'}->{'R2'}: $!";
	open($fh2, $r2OpenCmd) or die "Cannot open file $options{2}: $!";
}


#############################################################################
#read barcode data:
open(my $barFh, $options{'barcodes'}) or die "Cannot open file $options{'barcodes'}: $!";

while(<$barFh>){
	if(/^(\w+)\s+([ATGCN]+)\s*$/){
		my $name = $1;
		my $bar = $2;
		
		# print length($bar),"\t$barTrim\n";		
		if(length($bar) < $barTrim){
			print STDERR "Error: barcode length less than suggested barTrim option for barcode: ", $bar, "\n";			
			die;
		}
		
		## for 1-12 barcodes: use 7bp, for remaining barcodes: use 8bp
		$bar = substr(uc $bar, 0, $barTrim);
		
		#check if the same file name is used for different barcodes
		if(exists $files{$name}){
			print STDERR "Error: Two different barcodes are using same sample name.
			$files{$name} : $name
			$bar : $name\n";
			&error_cleanup();
			die;
		}
		
		$files{$name} = $bar;
		
		## check if a barcode is specified for two samples
		if(exists $barcodes{$bar}->{'name'}){
			print STDERR "Error: barcode clash. Two barcodes specified for same sample.
			$bar : $barcodes{$bar}->{'name'}
			$bar : $name\n";
			&error_cleanup();
			die;
		}
		

		$barcodes{$bar}->{'name'} = $options{'suffix'} eq '' ? $name : $name.'_'.$options{'suffix'};			#Sample name suffix
		$barcodes{$bar}->{'R1'} = $barcodes{$bar}->{'name'}.'_R1.fastq.gz';
		$barcodes{$bar}->{'R2'} = $barcodes{$bar}->{'name'}.'_R2.fastq.gz';
		
		#R1 file
		open($barcodes{$bar}->{'fhR1'}, "|-", "gzip > $barcodes{$bar}->{'R1'}") or die 'Cannot create file ',$barcodes{$bar}->{'R1'},": $!";
		
		if($isPaired){
			#R2 file
			open($barcodes{$bar}->{'fhR2'}, "|-", "gzip > $barcodes{$bar}->{'R2'}") or die "Cannot create file ",$barcodes{$bar}->{'R2'},": $!";
		}
		
		#this barcode length will be used to trim the sequence and qual line
		$barcodes{$bar}->{'barLen'} = length($bar);
		$barcodes{$bar}->{'count'} = 0;			#counter for each barcode
	}
	elsif(/^\s*$/){
		next;
	}
	else{
		&error_cleanup();
		print STDERR "Wrong format in barcode file at line: $_";
		die;
	}
}


#############################################################################
## demultiplexing

my ($p1, $p2, $outBar);

if($isPaired){
	#for paired end data
	while(1){
		if(eof($fh1) || eof($fh2)){
			if(eof($fh1) && eof($fh2)){
				last;
			}
			else{
				&error_cleanup();
				die "ERROR: Number of reads in R1 and R2 files are not equal\n";
			}
		}
		
		#read line1: headers
		$p1 = <$fh1>;
		$p2 = <$fh2>;
		
		#read line 2: sequence
		my $sq1 = <$fh1>;
		my $sq2 = <$fh2>;
		
		my $sq1Bar = substr($sq1, 0, $barTrim);
		my $sq2Bar = substr($sq2, 0, $barTrim);

		if(exists($barcodes{$sq1Bar})){
			$p1 .= $sq1;
			$p2 .= $sq2;
		}
		# elsif(exists($barcodes{$sq2Bar})){
			# $p1 .= $sq1;
			# $p2 .= $sq2;
		# }
		else{
			$sq1Bar = 'unknown';
			$p1 .= $sq1;
			$p2 .= $sq2;
		}
				
		#read line 3: +
		$p1 .= <$fh1>;
		$p2 .= <$fh2>;
		
		#read line 4: qual
		$p1 .= <$fh1>;
		$p2 .= <$fh2>;
		
		print {$barcodes{$sq1Bar}->{'fhR1'}} $p1;
		print {$barcodes{$sq1Bar}->{'fhR2'}} $p2;
		$barcodes{$sq1Bar}->{'count'}++;
	}
}
else{
	#for single end
	while(1){
		if(eof($fh1)){
			last;
		}
		
		#read line1: headers
		$p1 = <$fh1>;
		
		#read line 2: sequence
		my $sq1 = <$fh1>;

		my $sq1Bar = substr($sq1, 0, $barTrim);

		if(exists($barcodes{$sq1Bar})){
			$p1 .= $sq1;
		}
		else{
			$sq1Bar = 'unknown';
			$p1 .= $sq1;
		}
				
		#read line 3: +
		$p1 .= <$fh1>;
		
		#read line 4: qual
		$p1 .= <$fh1>;
		
		print {$barcodes{$sq1Bar}->{'fhR1'}} $p1;
		$barcodes{$sq1Bar}->{'count'}++;
					
	}
}

open(my $out, '>','demultiplex.stats') or die "Cannot create file reads.stats: $!";


foreach(sort{$barcodes{$a}->{'name'} cmp $barcodes{$b}->{'name'}}keys %barcodes){
	close($barcodes{$_}->{'fhR1'});
	if($isPaired){
		close($barcodes{$_}->{'fhR2'});
	}
	
	print $out $barcodes{$_}->{'name'}, "\t", $barcodes{$_}->{'count'}, "\n";
}

close($out);
close($barFh);



## cleanup the files in case of any error
sub error_cleanup{
	foreach(keys %barcodes){
		close($barcodes{$_}->{'fhR1'});
		unlink($barcodes{$_}->{'R1'});
		if($isPaired){
			close($barcodes{$_}->{'fhR2'});
			unlink($barcodes{$_}->{'R2'});
		}
	}
}














__END__


=head1 NAME


=head1 SYNOPSIS

perl demultiplex_chris_barcodes.pl --barcodes <Barcode file> --1 <Mate1> --2 <Mate2> --suffix <suffix to add>

Help Options:

	--help	Show this scripts help information.

=head1 DESCRIPTION

This script is used to demultiplex the libraries generated by Chris' multiplexing
protocol for ChIPseq samples.


=head1 OPTIONS

=over 30

=item B<--barcodes>

[STR] TAB separated Barcode file where first column is sample name and second column in barcode sequence

=item B<--barTrim>

[INT] Trim all barcodes to this length. Default: 8

=item B<--1>

[STR] Forward read file. COMMA separated multiple files can be provided.

=item B<--2>

[STR] Reverse read file (Optional). COMMA separated multiple files can be provided.

=item B<--suffix>

[STR] Sample name suffix to add (Optional)

=item B<--help>

Show this scripts help information.

=back


=cut

