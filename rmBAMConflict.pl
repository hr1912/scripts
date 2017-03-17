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
	my @line = split;
	next if ((!($line[1] & 0x04) && $line[2] eq "*") || (!($line[1] & 0x04) && $line[6] eq "*"));
	print OUT "$_";
}

close BAM; 
close OUT;
