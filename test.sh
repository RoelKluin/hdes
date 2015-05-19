
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
' | valgrind --leak-check=full ./uqct 2>&1 | ./fqless

#odd
make clean && make && zcat SRR077487_1.filt.fastq.gz | head -n10000 | valgrind --leak-check=full ./uqct 2>&1 | egrep -A3 \
"CAACATGGAGAAACCCC|GGGGTTTCTCCATGTTG" --color=always

make clean && make && zcat SRR077487_1.filt.fastq.gz | head -n1000000 | ./uqct 2>&1 | ./fqless
# search for CTCTGTGGTGTCTGATT

samtools view HG00096.chrom11.ILLUMINA.bwa.GBR.exome.20120522.bam "11:194121-194356" | perl -e '
while (<>) {
    my @L = (split /\t/)[1,9,2,10];
    if ($L[2] & 16) {
        $L[$_] = (scalar reverse $L[$_]) for (1, 3);
        $L[1] = join("", map { $_ =~ tr/ACGTacgt/TGCAtgca/; $_ } split(//, $L[1]));
    }
    $L[2] = "+";
    $L[0] = "@".$L[0];
    print join ("\n", @L)."\n";
}' | ./uqct | ./fqless


make clean && make && diff -u0 <(zgrep -o "^.* " head_100K_SRR077487_1.filt.fastq.gz | sort) <(./uqct head_100K_SRR077487_1.filt.fastq.gz | grep -o "^.* " | sort)


./uqct SRR081241_1.filt.fastq.gz -o bar.gz -f

diff -u0 <(zgrep -o "^.* " head_100K_SRR077487_1.filt.fastq.gz | sort) \
    <(./uqct head_100K_SRR077487_1.filt.fastq.gz | grep -o "^.* " | sort)

make clean && make
zgrep -E -n -A3 "^($(diff -u0 <(zgrep -o "^.* " head_100K_SRR077487_1.filt.fastq.gz | sort) \
    <(./uqct head_100K_SRR077487_1.filt.fastq.gz | grep -o "^.* " | sort) |
sed -n '/^-@/{s/^.//;H;}; ${x; s/\n/|/g;s/^|//p;}'))" head_100K_SRR077487_1.filt.fastq.gz


#./uqct head_100K_SRR077487_1.filt.fastq.gz | zgrep -E -n -A3 "^($(diff -u0 <(zgrep -o "^.* " head_100K_SRR077487_1.filt.fastq.gz | sort) \
#    <(./uqct head_100K_SRR077487_1.filt.fastq.gz | grep -o "^.* " | sort) |
#sed -n '/^-@/{s/^.//;H;}; ${x; s/\n/|/g;s/^|//p;}'))"



make clean && make && zcat /net/RTAdump/HiSeq2000/140716_M00872_0140_000000000-A9AAV_M153/Data/Intensities/BaseCalls/Unaligned/Project_MM153_Sander_Kelderman_2852/Sample_2852_29_CRC6_1_C5/2852_29_CRC6_1_C5_ATTGAGGA_L001_R1_001.fastq.gz | head -n10000 | valgrind ./uqct 200 2>&1 | ./fqless

make clean && make && ./uqct hg19.fa.gz

#zcat hg19_including_nonprimary.fa.gz |
#awk '/^>/ {FA="fa_parts/" substr($0,2,index($0," ")-1)}; {print | "gzip -c --fast - >> "FA; close(FA)}'

#for chr in `seq 1 22` X Y MT; do
#    zcat fa_parts/$chr
#done | gzip --fast > hg19.fa.gz

samtools faidx hg19.fa.gz
cut -f 1,2 hg19.fa.gz.fai | uniq | sed 's/^/chr/' > chrom.sizes

make clean && make && ./uqct hg19.fa.gz -l 51 2>&1 | tee uqct.err
zcat hg19_kcpos.wig.gz | ucsc/userApps/bin/wigToBigWig -clip stdin chrom.sizes /tmp/out.bw
~/dev/git/IGV/igv.sh /tmp/out.bw

zgrep -P "(0\s+){7}" /tmp/out.bw.gz -c
#203590084
echo $((1<<29))
#536870912
#=> 0.3792 empty


zcat /tmp/out.bw.gz | awk 'BEGIN {
 for (x = 0; x != 5; ++x) m[x] = 0;
} {
    if ($3 != "maxdbit") for (x = 0; x != 5; ++x) if ($(x+3) > m[x]) m[x]=$(x+3);
} END {
    for (x = 0; x != 5; ++x) print "column " x + 3 ":\t" m[x];
}' | tee maxes.txt

######################################
# smaller ref, note key is now smaller as well
#zcat hg19_including_nonprimary.fa.gz | sed '/^>[^GM]/Q' | gzip --fast > hg19_GL.fa.gz
make clean && DEFINES="-DKEY_LENGTH=11" make && ./uqct hg19_GL.fa.gz -l 51 2>&1 | tee uqct.err

./faless hg19_GL.fa.gz --extend 10



# test reading 1st stored file
(rm hg19_GL.{2b,nn,bd,ub,kc}.gz; make clean && DEFINES="-DKEY_LENGTH=11" make &&
    valgrind ./uqct hg19_GL.fa.gz -l 51 &&
    valgrind ./uqct hg19_GL.fa.gz -l 51)  2>&1 | tee hg19_GL_part1.err


rm hg19_GL.{2b,nn,bd,ub,kc}.gz; make clean && DEFINES="-DKEY_LENGTH=11" make &&
    ./uqct hg19_GL.fa.gz -l 51 && rm hg19_GL.ub.gz &&
    valgrind ./uqct hg19_GL.fa.gz -l 51 2>&1 | tee hg19_GL_part1.err

rm hg19_GL.x{1,2}.gz; make clean && DEFINES="-DKEY_LENGTH=11" make &&
 ./uqct hg19_GL.fa.gz -l 51 &> /dev/null && rm hg19_GL.x2.gz &&
 valgrind ./uqct hg19_GL.fa.gz -l 51 -b 64 2>&1 | less

rm hg19_GL.x{1,2}.gz; make clean && DEFINES="-DKEY_LENGTH=15" make &&
 ./uqct hg19_GL.fa.gz -l 51 &> /dev/null && rm hg19_GL.x2.gz &&
 valgrind ./uqct hg19_GL.fa.gz -l 51 -b 64 2>&1 | less

rm hg19_chr1.x1.gz; make clean && DEFINES="-DKEY_LENGTH=15" make && ./uqct hg19_chr1.fa.gz -l 51 2>&1 | tee hg19_chr1_uqct.err

bug bij 2e N stretch: niet gestopt voor stretch 

./faless hg19_GL.fa.gz --extend 10



# the real thing
#rm hg19_x1.gz;
make clean && DEFINES="-DKEY_LENGTH=15" make && ./uqct hg19.fa.gz -l 51 2>&1 | tee hg19_uqct.err


zcat hg19.fa.gz | sed '/^>3/Q' | gzip --fast > hg19_chr1_2.fa.gz

rm test.x*.gz; DEFINES="-DKEY_LENGTH=11" make && ./uqct test.fa.gz -l 51 2>&1 | tee test_uqct.err

make clean && DEFINES="-DKEY_LENGTH=15" make && valgrind ./uqct hg19_chr1_2.fa.gz -l 51 2>&1 | tee hg19_1_2uqct.err



