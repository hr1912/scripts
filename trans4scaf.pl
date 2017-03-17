#! /ngs/SD/ruanhang/ActivePerl-5.14/bin/perl -w
use strict;
use Getopt::Long;
use Cwd;
use Math::Random::OO::Normal;
use Data::Dumper;

=head
---------------------------------------------------------------
This script gets simulated contig sequences and generate associated files for iteration of soapdenovo scaffolding.

Author: Ruan Hang <ruanhang@novogene.cn>
----------------------------------------------------------------
=cut

if (@ARGV < 2) {
  die "Usage: trans4scaf.pl -k kmer_length -c haploid_coverage -g contig_file -o output_prefix\n";
}

my ($overlaplen,$h_cov) = (0,0);
my ($contig,$newcontig,$ctgindex,$upedge,$prefix) = ('','','','','');

GetOptions ('k=i' => \$overlaplen,
			'c=i' => \$h_cov,
			'g=s' => \$contig,
			'o=s' => \$prefix);

my $dir = getcwd;
`mkdir $dir/scaff`;
`touch $dir/scaff/$prefix.Arc`;
$newcontig = "$dir/scaff/$prefix.contig";
$ctgindex = "$dir/scaff/$prefix.ContigIndex";
$upedge = "$dir/scaff/$prefix.updated.edge";

open(CTG,$contig) || die "Incorrect file $contig, now exiting...\n";
open(NEW,">$newcontig") || die "Can not create file: ($!)";
open(IDX,">$ctgindex") || die "Can not create file: ($!)";
open(EDG,">$upedge") || die "Can not create file: ($!)";

my $timeS = time;
my %hash;

my $num_ctg = 0;

my ($id,$pos,$len,$cat,$seq) = (0,0,0,'','');
while(<CTG>){
	chomp;
	my $line = $_;
	$seq.= $line if (eof(CTG));
	if (/^>/ || eof(CTG)) {
		if ($seq ne '') {
			$hash{"$id:$pos"} = "$cat:$len:$seq";
			#($id,$pos,$len,$cat,$seq) = (0,0,0,'','');
			#$seq = '';
		}
		$seq = ''; 
		($id,$pos,$len,$cat) = (split /:/,$_)[0..3];
		$num_ctg ++;
		next;
	}
	$seq.= $line;
}

close CTG;
#print Dumper(\%hash);

my $time1 = time;
printf "Contig sequence file split and hashed~ Time consumed:%.2fmin\n", ($time1 - $timeS)/60;

print "Total number of contigs: $num_ctg\n";
printf IDX "Edge_num %d %d\n", $num_ctg*2, $num_ctg; 
print IDX "index\tlength\treverseComplment\n";
printf EDG "EDGEs %d\n", $num_ctg * 2; 

my $num_edge = $num_ctg*2;
`echo "VERTEX 100000 K $overlaplen\n\nEDGES $num_edge\n\nMaxReadLen 100 MinReadLen 0 MaxNameLen 256" > $dir/scaff/$prefix.PreGraphBasic`;

my $counter=1;
my $Normal_dis = Math::Random::OO::Normal->new($h_cov,20);
$Normal_dis->seed(time ^ $$);
foreach my $key (sort by_specific keys %hash){
	my @val = (split /:/, $hash{$key})[0,1,2];
	my $oneCvg=0;
	my $category = 0;
	if ($val[0] eq 'Hete') {
		$oneCvg = sprintf("%.1f",$Normal_dis->next());
		$category = 1;
	} elsif ($val[0] eq 'Homo'){
		$oneCvg = sprintf("%.1f",$Normal_dis->next()*2);
		$category = 2;
	} else {
		$oneCvg = sprintf("%.1f",$Normal_dis->next()*(rand(3)+3));
		$category = 4;
	}
	$oneCvg = 0 if ($oneCvg < 0);
	printf NEW ">%d length %d cvg_%.1f_tip_0_cat_%d\n",$counter,$val[1], $oneCvg,$category;
	print NEW "$val[2]\n";
	printf IDX "%d\t%d\t1\n",$counter,$val[1];
	printf EDG ">length %d,1,%d 0 0\n",$val[1]-$overlaplen,$oneCvg*10;
	printf EDG ">length %d,-1,%d 0 0\n",$val[1]-$overlaplen,$oneCvg*10;
	$counter += 2;
}

close NEW;
close IDX;
close EDG;

my $timeE = time;
printf "Files for scaff generated~ Time consumed:%.2fmin\n", ($timeE - $time1)/60;
printf "All done~ Total time consumed:%.2fmin\n", ($timeE - $timeS)/60;

sub by_specific {
	(split /:/,$hash{$a})[1] <=> (split /:/,$hash{$b})[1]
		or
	(split /:/,$hash{$a})[0] cmp (split /:/,$hash{$b})[0]
		or
	(split /:/,$hash{$a})[2] cmp (split /:/,$hash{$b})[2]; 
}

