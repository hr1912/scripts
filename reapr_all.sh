PATH=/ngs/SD/ruanhang/ActivePerl-5.14/site/bin:/ngs/SD/ruanhang/ActivePerl-5.14/bin:/ngs/self-software/gcc-4.7.0/bin:$PATH
export PATH
LD_LIBRARY_PATH=/ngs/SD/ruanhang/zlib-1.2.7/lib:/ngs/self-software/gcc-4.7.0/lib64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
C_INCLUDE_PATH=/ngs/SD/ruanhang/zlib-1.2.7/include:$C_INCLUDE_PATH
export C_INCLUDE_PATH

export SMALT="/database/ruanhang/bin/smalt"
export BWA="/ngs/self-software/bin/bwa"
export SAM="/ngs/self-software/bin/samtools"
export REAPR="/database/ruanhang/bin/reapr" 

$REAPR facheck Hum_chr14.k35.scafSeq assembled
sed -i '/^$/d' assembled.fa
#$SMALT index -k 13 -s 6 smalt_index SOAPassembled.fa
#$SMALT map -n 12 -f sam -o frag.sam -x smalt_index ../data/frag_1.fastq ../data/frag_2.fastq 
#$SMALT map -n 12 -f sam -o short_jump.sam -x smalt_index ../data/shortjump_1.fastq ../data/shortjump_2.fastq

$BWA index -p bwa_index assembled.fa
#$BWA aln -t 12 bwa_index frag_1.fastq > frag_1.sai &
#$BWA aln -t 12 bwa_index frag_2.fastq > frag_2.sai
#$BWA sampe -a 250 -P bwa_index frag_1.sai frag_2.sai frag_1.fastq frag_2.fastq | $SAM view -bS - | $SAM sort -m 15000000000 - frag.sorted
 
#$SAM sort -m 15000000000 frag.bam frag.sorted
#$SAM rmdup frag.sorted.bam frag.sorted.rmdup.bam | /database/ruanhang/bin/processBam4Reapr.pl - 
#/database/ruanhang/bin/processBam4Reapr.pl frag.sorted.rmdup.bam frag.final.bam
#$SAM index frag.final.bam

$BWA aln -t 12 bwa_index shortjump_1.fastq > shortjump_1.sai &
$BWA aln -t 12 bwa_index shortjump_2.fastq > shortjump_2.sai
$BWA sampe -a 3500 -P bwa_index shortjump_1.sai shortjump_2.sai shortjump_1.fastq shortjump_2.fastq | $SAM view -bS - | $SAM sort -m 15000000000 - frag.sorted

#$SAM sort -m 15000000000 short_jump.bam short_jump.sorted
$SAM rmdup short_jump.sorted.bam short_jump.sorted.rmdup.bam
/database/ruanhang/bin/processBam4Reapr.pl short_jump.sorted.rmdup.bam short_jump.final.bam
$SAM index short_jump.final.bam

#$REAPR perfectfrombam frag.final.bam perfect 100 260 3 4 76
$REAPR pipeline assembled.fa short_jump.final.bam reapr_out perfect
