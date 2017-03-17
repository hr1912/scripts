#! /ngs/SD/ruanhang/ActivePerl-5.14/bin/perl -w
use strict;
use 5.010001;
use PerlIO::gzip;
=head
---------------------------------------------------------------------

Author: Ruan Hang <ruanhang@novogene.cn>
---------------------------------------------------------------------
=cut

my $read_fp = shift or die "dumpReadsN.pl <reads.fq.gz> <reads.noN.fq.gz>\n";
my $new_fp = shift or die "dumpReadsN.pl <reads.fq.gz> <reads.noN.fq.gz>\n";

open (RDS,"<:gzip", "$read_fp") || die "$!";
open (NEW,">:gzip", "$new_fp") || die "$!";

my ($tol_rd_count,$rd_n_count,$tol_n_count) = (0,0,0);

while (<RDS>) {
	if (/^@/) {
		my $line_1 = $_;
		my $line_2 = <RDS>;
		my $line_3 = <RDS>;
		my $line_4 = <RDS>;
		my $n_count = 0;
		$n_count += ($line_2 =~ tr/N/N/);
		$tol_rd_count ++;
		if ($n_count==0) {
			print NEW $line_1.$line_2.$line_3.$line_4;
		} else {
			$rd_n_count ++;
			$tol_n_count += $n_count;
		}
	}
}

close RDS;
close NEW;
print "$tol_rd_count,$rd_n_count,$tol_n_count\n";

