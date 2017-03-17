#! /ngs/SD/ruanhang/ActivePerl-5.14/bin/perl -w
use strict;
use Getopt::Long;
use Cwd;
use Data::Dumper;

=head
---------------------------------------------------------------
This script splits scaffold sequences into contig sequences and generate associated files for iteration of soapdenovo scaffolding.

Author: Ruan Hang <ruanhang@novogene.cn>
----------------------------------------------------------------
=cut

if (@ARGV < 2) {
  die "Usage: scaf4ite.pl <-r> -k kmer_length -g soapdenovo_prefix \n";
}
my $overlaplen = 0;
my $repeat_exclude;
my ($prefix,$scafseq,$gapseq,$ctgalonggap,$scafcov,$contig,$newcontig,$ctgindex,$upedge) = ('','','','','','','','','');

GetOptions ('k=i' => \$overlaplen,
			'g=s' => \$prefix,
			'r' => \$repeat_exclude);

my $dir = getcwd;
`mkdir $dir/ite`;
`touch $dir/ite/$prefix.Arc`;
if (-e "$dir/$prefix.preGraphBasic"){
	`cp $dir/$prefix.preGraphBasic $dir/ite`;
}
#system "touch $dir\/ite\/$contig.indexTrack";
$scafseq = "$dir/$prefix.scafSeq";
$gapseq = "$dir/$prefix.gapSeq";
$scafcov = "$dir/$prefix.scaf_cov";
$ctgalonggap = "$dir/$prefix.ctgAlongGap";
$newcontig = "$dir/ite/$prefix.contig";
$ctgindex = "$dir/ite/$prefix.ContigIndex";
$upedge = "$dir/ite/$prefix.updated.edge";
$contig = "$dir/$prefix.scaf_cut";

open(COV,$scafcov) || die "Incorrect file $scafcov, now exiting...\n";
open(CTG,$scafseq) || die "Incorrect file $scafseq, now exiting...\n";
if ($repeat_exclude) {
	open(GAP,$gapseq) || die "Incorrect file $gapseq, now exiting...\n";
	open(CAG,$ctgalonggap) || die "Incorrect file $ctgalonggap, now exiting..\n";
}
open(CUT,">$contig") || die "Can not create file: ($!)";
open(NEW,">$newcontig") || die "Can not create file: ($!)";
open(IDX,">$ctgindex") || die "Can not create file: ($!)";
open(EDG,">$upedge") || die "Can not create file: ($!)";

my %hash_c;
my ($timeS,$time1,$time2,$time3,$time4,$timeE);

$timeS = time;

my $id = '';
my ($cov,$cov_sum,$cov_n,$cov_avg) = (0,0,0,0);

while(<COV>){
	chomp;
	s/^\s+//g;
	if (/^>(\S+)/){
		$id = $1;
		chomp ($cov = <COV>);
		$cov =~ s/^\s+//g;
		push @{$hash_c{$id}}, $cov;
		$cov_n ++;
		$cov_sum += $cov;
	} else {
		push @{$hash_c{$id}}, $_;
		$cov_n ++;
		$cov_sum += $cov;
	}
}

#print Dumper(\%hash);
close COV;
$cov_avg = sprintf("%.1f",$cov_sum / $cov_n);

$time1 = time;
printf "Coverage Info attached~ Average contig coverage is %.1f\nTime consumed:%.2fmin\n", $cov_sum / $cov_n, ($time1 - $timeS)/60;

my %hash;

my ($seq, $name)=('', '');
my $len = 0;
my ($num_ctg,$num_ctg_tmp) = (0,0);

while(<CTG>){
  chomp;
  my $line = $_;
  $seq.= uc($line) if(eof(CTG));
  if (/^>(\S+)/ || eof(CTG)){
    if($seq ne ''){
      my @seqgaps = split(/[N]{1,}/, $seq);
      if($#seqgaps > 0){
        my $ctgcount=0;
        foreach my $ctgseq (@seqgaps){
		  #$num_ctg ++;
          #$ctgcount++;
		  $len = length($ctgseq);
		  next if ($len < $overlaplen);
		  $num_ctg ++;
		  $ctgcount ++;	
		  #print "$name\t$ctgcount\n";
		  #print "$hash{$name}\n";
          print CUT ">$name.$ctgcount\t$hash_c{$name}->[$ctgcount-1]\n$ctgseq\n";
		  $hash{"$name.$ctgcount"} = "$len:$hash_c{$name}->[$ctgcount-1]:$ctgseq";
        }
      }else{
		$num_ctg ++;
		#print "$name\n";
		$len = length($seq);
        print CUT ">$name\t$hash_c{$name}->[0]\n$seq\n";
		$hash{"$name"} = "$len:$hash_c{$name}->[0]:$seq";
      }
    }
    $seq='';
    #my @temp = split;
    #$name = $temp[0];
	if (/^>(\S+)/) {
		$name = $1;
	}
  }else{
    $seq.= uc($line);
  }
}

close CTG;
close CUT;

$time2 = time;
printf "Scaffold sequence file split and hashed~ %d contigs added.\nTime consumed:%.2fmin\n", $num_ctg,($time2 - $time1)/60;

if ($repeat_exclude) {
	$num_ctg_tmp = $num_ctg;

	$len = 0;
	($seq,$name) = ('','');

	while (<GAP>) {
		chomp;
		my $line = $_;
		$seq.= uc($line) if (eof(GAP));
		if (/^>S/ || eof(GAP)) {
			#print "*\t";
			if ($seq ne '') {
				$name = $1 if (/^>(S\d+?)/); 
				#print "*\t";
				$num_ctg ++;
				$len = length($seq);
				if ($len >= $overlaplen){
					#print "*\t";
					$hash{"$name"} = "$len:$cov_avg:$seq";
				}
			}
			$seq = '';
		} else {
			$seq .= uc($line);
		}	
	}

	close GAP;
	$time3 = time;
	printf "Gap sequence file hashed~ %d contigs added.\nTime consumed:%.2fmin\n", $num_ctg - $num_ctg_tmp, ($time3 - $time2)/60;

	$num_ctg_tmp = $num_ctg;

	$len = 0;
	($seq,$name) = ('','');

	while (<CAG>) {
		chomp;
		my $line = $_;
		$seq.= uc($line) if (eof(GAP));
		if (/^>(\d+)/ || eof(GAP)){
			if ($seq ne ''){
				$name = $1 if (/^>(\d+)/);
				$num_ctg ++;
				$len = length($seq);
				if ($len >= $overlaplen){
					$hash{"$name"} = "$len:$cov_avg:$seq";
				}	
			}
			$seq = '';
		} else {
			$seq .= uc($line);
		}
	}
	close CAG;
	$time4 = time;
	printf "ctgAlongGap sequence file hashed ~ %d contigs added.\nTime consumed:%2.fmin\n", $num_ctg - $num_ctg_tmp,($time4 - $time3)/60;

}

print "Total number of contigs: $num_ctg\n";
printf IDX "Edge_num %d %d\n", $num_ctg*2, $num_ctg; 
print IDX "index\tlength\treverseComplment\n";
printf EDG "EDGEs %d\n", $num_ctg * 2; 

my $counter=1;

foreach my $key (sort by_specific keys %hash) {
	my $flag=0;
	$flag = 1 unless ($key =~ /scaffold/);
	my @val = (split /:/, $hash{$key})[0,1,2];
	print NEW ">$counter length $val[0] cvg_$val[1]_tip_$flag\n";

	$val[2] =~ s/\s+//g;
		
	for (my $i = 0; $i < length($val[2]); $i += 100) {
		my $string = substr($val[2],$i,100);
		print NEW $string,"\n";
	} 
	
	my $flag_palindrome = 1;
	
	if ($val[0] % 2 == 0) {
		$flag_palindrome = &isPalindrome($val[2]);
	}
	
	print IDX "$counter\t$val[0]\t$flag_palindrome\n";
	if ($flag_palindrome == 0) {
		printf EDG ">length %d,0,%d 0 0\n",$val[0]-$overlaplen,$val[1]*10;
		$counter ++;
	} elsif ($flag_palindrome == 1) {
		printf EDG ">length %d,1,%d 0 0\n",$val[0]-$overlaplen,$val[1]*10;
		printf EDG ">length %d,-1,%d 0 0\n",$val[0]-$overlaplen,$val[1]*10;
		$counter += 2;
	}
}

close NEW;
close IDX;
close EDG;

$timeE = time;
if ($repeat_exclude) {
	printf "Files for iteration generated~ Time consumed:%.2fmin\n", ($timeE - $time4)/60;
} else {
	printf "Files for iteration generated~ Time consumed:%.2fmin\n", ($timeE - $time2)/60;
}
printf "All done~ Total time consumed:%.2fmin\n", ($timeE - $timeS)/60;

sub by_specific {
	(split /:/,$hash{$a})[0] <=> (split /:/,$hash{$b})[0]
		or
	$a cmp $b;
}

sub isPalindrome {
	my $str = $_[0];
	my @seq = split(//,$str);
	my $i = 0;
	my ($head, $tail) = ('', '');
	my $len = scalar(@seq);
	print $len if ($len % 2 != 0) ;
	while (1) {
		$i ++;
		$head = uc(shift(@seq));
		$tail = &ComplementBase(pop(@seq));
		return 1 if ($head ne $tail);
		return 0 if ($i >= $len/2);
	}
}

sub ComplementBase {
	my $char = uc($_[0]);
	$char =~ tr/ACTG/TGAC/;
	return $char;
}
