#! /ngs/SD/ruanhang/ActivePerl-5.14/bin/perl -w
use strict;
use Cwd;
use Getopt::Long;
#use Data::Dumper;
=head
---------------------------------------------------------------------

Author: Ruan Hang <ruanhang@novogene.cn>
---------------------------------------------------------------------
=cut

if (@ARGV <2) {
	die "Usage: prepare4scaff.pl -k kmer_length -c contig_original -g soapdenovo_prefix \n";
}

my ($overlaplen,$num_ctg) = (0,0);
my ($prefix, $contig_fp,$ctgindex_fp,$upedge_fp,$newctg_fp) = ('','','','','');

GetOptions ('k=i' => \$overlaplen,
			'c=s' => \$contig_fp,
			'g=s' => \$prefix);

my $dir = getcwd;
`mkdir $dir/files4scaff`;
`touch $dir/files4scaff/$prefix.Arc`;

=head
if (-e "$dir/$prefix.preGraphBasic") {
	`cp $dir/$prefix.preGraphBasic $dir/scaff`;
} 
=cut

$ctgindex_fp = "$dir/files4scaff/$prefix.ContigIndex";
$upedge_fp = "$dir/files4scaff/$prefix.updated.edge";
$newctg_fp = "$dir/files4scaff/$prefix.contig";

open(CTG,$contig_fp) || die "$!";
open(IDX,">$ctgindex_fp") || die "$!";
open(EDG,">$upedge_fp") || die "$!";
open(NEW,">$newctg_fp") || die "$!";

my %ctgseq;

my ($seq, $lastname, $name) = ('','','');

while(<CTG>) {
	chomp;
	#$seq .= $_ if (eof(CTG));
	if (/^>(\S+)/ || eof(CTG)) {
		$name = $1;
		$ctgseq{$lastname} = $seq if ($lastname);
		$lastname = $name;  
		$num_ctg ++;
		$seq = '';
	} else {
		$seq .= $_;
	}
}

#print Dumper(\%ctgseq);

close CTG;

open(CTG,$contig_fp) || die "$!";

printf IDX "Edge_num %d %d\n", $num_ctg*2, $num_ctg;
print IDX "index\tlength\treverseComplment\n";
printf EDG "EDGEs %d\n", $num_ctg * 2;

my $num_edge = $num_ctg*2;
`echo "VERTEX 100000 K $overlaplen\n\nEDGES $num_edge\n\nMaxReadLen 100 MinReadLen 0 MaxNameLen 256" > $dir/files4scaff/$prefix.preGraphBasic`;

my $counter = 1;
while(<CTG>) {
	chomp;
	if (/^>/) {
		my @line = split;
		my $id = $line[0];
		$id =~ s/^>//g;
		my $len = $line[2];
		my $cvg;
		if ($line[3] =~ /cvg_(\S+)_tip/){
			$cvg = $1;
		}
		my $flag = 1;

		if ($len % 2 == 0) {
			#print "$id\t\n"; 
			$flag = &CtgDiffFromBal($id);
		}
		
		#print "$ctgseq{$id}\n" if (!$flag);
		print NEW ">$counter $line[1] $line[2] $line[3]\n";		
		print IDX "$counter\t$len\t$flag\n";
		if ($flag == 1) {
			printf EDG ">length %d,1,%d 0 0\n",$len-$overlaplen,$cvg*10;
			printf EDG ">length %d,-1,%d 0 0\n",$len-$overlaplen,$cvg*10;
			$counter += 2;	
		} elsif ($flag == 0) {
			printf EDG ">length %d,0,%d 0 0\n",$len-$overlaplen,$cvg*10;
			$counter ++;	
		}
	} else {
		print NEW "$_\n";
	} 
}

close CTG;
close IDX;
close EDG;
close NEW;

sub CtgDiffFromBal {
	my $id = $_[0];
	#print "$id\n";
	if (!exists $ctgseq{$id}) {
		print STDERR "Contig ID do not exists\n";
		return -1;
	} 
	my @seq = split(//,$ctgseq{$id});
	#print "$ctgseq{$id}\n";
	my $i = 0;
	my ($head, $tail) = ('', '');
	my $len = scalar(@seq);
	while (1) {
		$i ++;	
		$head = uc(shift(@seq)) ;
		$tail = &ComplementBase(pop(@seq));
		return 1 if ($head ne $tail);
		return 0 if ($i > $len/2) ;
	} 
}

sub ComplementBase {
	my $char = uc($_[0]);
	$char =~ tr/ACTG/TGAC/;
	return $char;
}
