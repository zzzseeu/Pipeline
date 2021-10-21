#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw[min max];
use File::Basename;

my $usage = 'This script extract the top 10 frequent sequences in the begining of reads in the fastq file.
These frequent sequences are then crossmatched with the barcode list to see if any sequence
matches better to any barcode other than intended barcode.

USAGE: perl freq_barcode_crossmatch.pl <input_file>
<input_file> should be a TAB delimited file and should have three columns: <sample_id>\t<fastq_file_path>\t<barcode>
';



if(!defined($ARGV[0])){
	print STDERR $usage;
	die "\nERROR: Input file not provided\n";
}


my %barcodes = ();

foreach(<DATA>){
	chomp $_;
	my ($id, $seq) = split(/\t/, $_);
	$barcodes{$seq} = $id;
}



open(my $fh, $ARGV[0]) or die "File not found: $_[0] ", $!;

print join("\t", qw(sample barcodeId barcode totalReads kmerId kmer kmerCount fraction kmerMismatches crossmapBarId crossmapBar match crossmap)),"\n";

while(my $line = <$fh>){
	chomp $line;
	next if($line=~m/^\s*$/);
	
	my ($sp, $bc, $fq) = split(/\t/, $line);
	my $bcLen = length($bc);

	# print "$fq\t$bc\t$bcLen\n";
	my $i=3;
	my %freq=();
	
	## frequency of first n base pairs
	open(my $fastq, "gzip -dc $fq |") or die;	
	while(<$fastq>){
		
		# print $_;
		# last;

		if($i % 4 == 0){
			$freq{substr($_, 0, $bcLen)}++;
		}
		
		$i++;		
	}
	
	my $readCount = ($i - 3) / 4;
	my $count = 0;
	
	foreach my $ad(sort{$freq{$b} <=> $freq{$a}}keys %freq){
		$count++;

		my $adMis = &index_match($ad, $bc);
		my $matchPattern = ($ad ^ $bc) | ( "\x30" x max(length($ad), $bcLen) );
		$matchPattern =~ s/[^0]/*/g;
		print join("\t", $sp, $barcodes{$bc}, $bc, $readCount, "frag$count", $ad, $freq{$ad}, sprintf("%.3f", $freq{$ad}/$readCount), $adMis, $barcodes{$bc}, $bc, $matchPattern, "original"),"\n";
		
		foreach my $bar(sort keys %barcodes){
			my $barMis = &index_match($ad, $bar);
			if($barMis < $adMis && $bar ne $bc){
				my $matchPattern = ( $ad ^ $bar ) | ( "\x30" x max(length($ad), length($bar)) );
				$matchPattern =~ s/[^0]/*/g;
				print join("\t", $sp, $barcodes{$bc}, $bc, $readCount, "frag$count", $ad, $freq{$ad}, sprintf("%.3f", $freq{$ad}/$readCount), $barMis, $barcodes{$bar}, $bar, $matchPattern, "new"),"\n";
			}
		}
		
		
		if($count >= 10){
			last;
		}
	}

	close($fastq);
	
	# print "\n";
}

close($fh);










# # sample code for sequence mismatch count
# my $s1 = 'ATGCCATGC';
# my $s2 = 'ATGCAAT';

# my $x = $s1 ^ $s2;
# print "|",$x,"|\n";
# print $x | ("\x30" x length($s2)),"::\n";
# $x =~ tr/\0//d;
# print "|",$x,"|\n";
# print length($x),"**\n";
# my $y = $s1 =~ tr/AT//rd;
# print "||",$y,"||\n";

# print length($x),"\n";

sub index_match{
	my ($s1, $s2) = ('', '');
	
	## ensure that shorter sequence is $s1 and longer is $s2
	if(length($_[0]) < length($_[1])){
		$s1 = $_[0];
		$s2 = $_[1];
	}
	elsif(length($_[0]) > length($_[1])){
		$s1 = $_[1];
		$s2 = $_[0];
	}
	else{
		$s1 = $_[0];
		$s2 = $_[1];
	}
	
	
	my $mismatches = length($s1) - ( ( $s1 ^ $s2 ) =~ tr/\0// );
	
	return $mismatches;
}



## inhouse barcodes
__DATA__
bar01	ATCACGT
bar02	CGATGTT
bar03	TTAGGCT
bar04	TGACCAT
bar05	ACAGTGT
bar06	GCCAATT
bar07	CAGATCT
bar08	ACTTGAT
bar09	GATCAGT
bar10	TAGCTTT
bar11	GGCTACT
bar12	CTTGTAT
bar13	ATATAGGAT
bar14	AACCGTGTT
bar15	AGGTCAGTT
bar16	CTCTGTCTT
bar17	CCATACACT
bar18	CGCATTAAT
bar19	GTCTACATT
bar20	GAGTTAACT
bar21	GCAGCCTCT
bar22	TCGCGTACT
bar23	TATACCGTT
bar24	TGCGGTTAT
bar25	AACACCTACT
bar26	CCTTTACAGT
bar27	GGTCCTTGAT
bar28	TTGAGTGTT
bar29	ACTAACTGCT
bar30	CAGGAGGCGT
bar31	GTTGTCCCAT
bar32	TGACGCATT
bar33	ATCGCCAGCT
bar34	CATTCCAAGT
bar35	GCAAGTAGAT
bar36	TGATCCGAT
bar37	ACGTAGCTCT
bar38	CGAACTGTGT
bar39	TAGCTAGTAT
bar40	GTGGGATAT
bar41	ATCCTATTCT
bar42	CGGACGTGGT
bar43	GCGTTTCGAT
bar44	TATCTCCGT
bar45	ACAGTGCACT
bar46	CACAGTTGGT
bar47	GTGACTACAT
bar48	TGAGAGTGT
bar49	AATGCTGACT
bar50	CCGTCTGAGT
bar51	GGCAGACGAT
bar52	TTCTGATGT
bar53	AGTAGTGGCT
bar54	CTAGTCATGT
bar55	GACACTCTAT
bar56	TCATTAGGT
bar57	TCCAGCCTCT
bar58	CTAGATTCGT
bar59	GAACGCTGAT
bar60	AGAACACCT