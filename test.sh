
#make && zgrep -E -o "^[ACTGN]{6}" /home/roel/dev/rproject/rmap/hg19.fa.gz | head -n 800000 | valgrind ./uqct 2>&1 | less

#time zgrep -E -o "^[ACTGN]{6}" /home/roel/dev/rproject/rmap/hg19.fa.gz | head -n 800000 | ./uqct


make && zcat /home/roel/dev/git/bwa/my/bwa/test/ce_med.fq.gz | valgrind ./uqct 1>/dev/null 2> valg.err

make && zcat /home/roel/dev/git/bwa/my/bwa/test/ce_med.fq.gz | ./uqct | less

# test that entries in fastq remain the same
f=SRR081241_1.filt.fastq.gz

cmp <(zcat "$f" | head -n 400000 | sed -n '/[[:space:]]/s/[[:space:]].*$//;N;h;n;x;N;y/\n/\t/;p' | sort) \
    <(zcat "$f" | head -n 400000 | ./uqct | sed -n '/[[:space:]]/s/[[:space:]].*$//;N;h;n;x;N;y/\n/\t/;p' | sort) &&
        echo "$f raw and uqct output have the same fastq entries"

#should result in same deviant:
make clean && make && zcat SRR077487_1.filt.fastq.gz | head -n4 | perl -e '
my @e = (<>, <>, <>, <>);
print join("", @e);
print join("", $e[0], substr($e[1],1,-1)."A\n", $e[2], substr($e[3],1,-1)."H\n");
$e[$_] = (scalar reverse substr($e[$_], 0, -1))."\n" for (1, 3);
$e[1] = join("", map { $_ =~ tr/ACGTacgt/TGCAtgca/; $_ } split(//, $e[1]));
print join("", @e);
print join("", $e[0], substr($e[1],1,-1)."A\n", $e[2], substr($e[3],1,-1)."H\n");
' | valgrind --leak-check=full ./uqct 100 2>&1 | ./fqless

#odd
make clean && make && zcat SRR077487_1.filt.fastq.gz | head -n10000 | valgrind --leak-check=full ./uqct 2>&1 | egrep -A3 \
"CAACATGGAGAAACCCC|GGGGTTTCTCCATGTTG" --color=always

make clean && make && zcat SRR077487_1.filt.fastq.gz | head -n1000000 | ./uqct 100 2>&1 | ./fqless
# search for CTCTGTGGTGTCTGATT

samtools view HG00096.chrom11.ILLUMINA.bwa.GBR.exome.20120522.bam "11:194121-194356" | perl -e '
while (<>) {
    my @L (split /\t/)[1,9,2,10];
    if ($L[2] & 16) {
        $L[$_] = (scalar reverse $L[$_]) for (1, 3);
        $L[1] = join("", map { $_ =~ tr/ACGTacgt/TGCAtgca/; $_ } split(//, $L[1]));
    }
    $L[2] = "+";
    $L[0] = ">".$L[0];
    join ("\n", @L)."\n";
}'

