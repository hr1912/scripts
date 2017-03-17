#! /ngs/SD/ruanhang/ActivePerl-5.14/bin/perl -w
use strict;
use 5.010001;

=head
---------------------------------------------------------------------

Author: Ruan Hang <ruanhang@novogene.cn>
---------------------------------------------------------------------
=cut

my $infile = shift or die $!;
my $outfile = shift or die $!;
open BAM , "/ngs/self-software/bin/samtools view -h $infile |" or die $!;
open OUT , "| /ngs/self-software/bin/samtools view -bS - >$outfile" or die $!;

while (<BAM>) {
	if (/^@/) {
		print OUT "$_";
		next;
	}		
	chomp;
	my @line = split;
	next if (@line < 12);
	next if ((!($line[1] & 0x04) && $line[2] eq "*") || (!($line[1] & 0x04) && $line[6] eq "*"));
	my $cigar = $line[5];
	my $md = '';
	if ($line[$#line] =~ /MD/) {
		$md = (split /:/, $line[$#line])[2];
	} elsif ($line[$#line-1] =~ /MD/) {
		$md = (split /:/, $line[$#line-1])[2];
	} else {
		print STDERR "err!\n";
	}

	my ($match_plus_mis,$ins,$del,$md_del) = (0,0,0,0); 
	my ($match,$mismatch,$gap_opn,$gap_ext) = (0,0,0,0);
	
	while ($cigar =~ m/(\d+)M/g) {
		 $match_plus_mis += $1;
	}
	
	while ($cigar =~ m/(\d+)I/g) {
		$ins += $1;
		$gap_opn ++;
		$gap_ext += ($1-1); 
	}
	
	while ($cigar =~ m/(\d+)D/g) {
		$del += $1;
		$gap_opn ++;
		$gap_ext += ($1-1);
	}

	$mismatch += ($md =~ tr/ATCG/ATCG/);

	if ($del) {
		while ($md =~ m/\^([ACGT]+)/g) {
			$md_del += length($1) ;
		}
		$mismatch -= $md_del;
	}
	
	$match = $match_plus_mis - $mismatch;

	my $score = $match + $mismatch * (-2) + $gap_opn * (-4) + $gap_ext * (-3);
	$score = 0 if ($score < 0); 	
	print OUT "$_\tAS:i:$score\n";
}

close BAM; 
close OUT;
