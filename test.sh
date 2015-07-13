
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
./tst hg19_GL.fa.gz -v -c -m 11

./tst hg19_GL.fa.gz -p -m 11

./tst hg19_chr1_2.fa.gz -c -m 15
./tst hg19_chr1_2.fa.gz -p -m 15

./tst hg19.fa.gz -p -m 16



(rm hg19_GL.{2b,nn,bd,ub,kc}.gz; make clean && DEFINES="-DKEY_LENGTH=9" make &&
    valgrind ./uqct hg19_GL.fa.gz -l 51; rm hg19_GL.ub.gz; 
    valgrind ./uqct hg19_GL.fa.gz -l 51)  2>&1 | tee hg19_GL_part1.err | less

rm hg19_GL.ub.gz; .make clean && DEFINES="-DKEY_LENGTH=9" make && valgrind ./uqct hg19_GL.fa.gz -l 51 2>&1 | less

rm hg19_GL.ub.gz && make clean && DEFINES="-DKEY_LENGTH=11" make &&
valgrind ./uqct hg19_GL.fa.gz -l 51 2>&1 | tee hg19_GL_part1.err


rm hg19_GL.{2b,nn,bd,ub,kc}.gz; make clean && DEFINES="-DKEY_LENGTH=15" make &&
    ./uqct hg19_GL.fa.gz -l 51 && rm hg19_GL.ub.gz &&
    valgrind ./uqct hg19_GL.fa.gz -l 51 2>&1 | tee hg19_GL_part1.err

rm hg19_chr1.x1.gz; make clean && DEFINES="-DKEY_LENGTH=15" make &&
    ./uqct hg19_chr1.fa.gz -l 51 2>&1 | tee hg19_chr1_uqct.err

bug bij 2e N stretch: niet gestopt voor stretch 

./faless hg19_GL.fa.gz --extend 10



# the real thing
#rm hg19.{2b,nn,bd,ub,kc}.gz;
make clean && DEFINES="-DKEY_LENGTH=16" make && ./uqct hg19.fa.gz -l 51 2>&1 | tee hg19_uqct.err


zcat hg19.fa.gz | sed '/^>3/Q' | gzip --fast -c > hg19_chr1_2.fa.gz
make clean && DEFINES="-DKEY_LENGTH=15" make &&
valgrind ./uqct hg19_chr1_2.fa.gz -l 51 2>&1 | tee hg19_1_2uqct.err


rm test.x*.gz; DEFINES="-DKEY_LENGTH=11" make && ./uqct test.fa.gz -l 51 2>&1 | tee test_uqct.err

mkdir fakeq
for ref in hg19_GL.fa.gz; do
  for simerr in "" "-SE true"; do
    out="fakeq/${ref%.fa.gz}${simerr// /_}"
    echo $out
    /opt/java/jre1.8.0_45/bin/java \
-jar  /home/roel/dev/git/hdes/external/ArtificialFastqGenerator.jar \
-R <(zcat hg19_GL.fa.gz) -O $out -RL 51 -RCNF 1 -SE true \
-S $(zcat hg19_GL.fa.gz | head -n 1 | sed -r 's/^(>[^ \t]+)([ \t].*)?$/\1/') 2> /dev/null
  done
done


make clean && DEFINES="-DKEY_LENGTH=11" make &&
valgrind ./uqct fakeq/hg19_GL.1.fastq.gz hg19_GL.2b.gz -l 51 2>&1 | tee hg19_GL.samerr

valgrind ./uqct fakeq/hg19_GL.1.fastq.gz hg19_GL.2b.gz

#qualities are not 
#####################
#rm hg19_GL.{2b,nn,bd,ub,kc,uq}.gz;
make clean && DEFINES="-DKEY_LENGTH=11" make &&
valgrind ./uqct hg19_GL.fa.gz -l 51 2>&1 | tee hg19_GL_uqct.err
####################

#with debug 7
cd /home/roel/dev/git/hdes &&
make clean &&
DEFINES="-DKEY_LENGTH=11" make &&
cd /home/roel/dev/git/hdes/bwatest &&
valgrind ../uqct ../fakeq/hg19_GL.1.fastq.gz ../hg19_GL.2b.gz -l 51 2>&1 | less

cd /home/roel/dev/git/hdes
make clean
DEFINES="-DKEY_LENGTH=11" make
cd /home/roel/dev/git/hdes/bwatest

mkdir /home/roel/dev/git/hdes/bwatest; cd !$
ln -s /home/roel/dev/git/hdes/hg19_GL.fa.gz
bwa=/net/NGSanalysis/apps/bwa/bwa-0.7.12/bwa
samtools=/net/NGSanalysis/apps/samtools/samtools-0.1.19/samtools
$bwa index hg19_GL.fa.gz
$samtools faidx hg19_GL.fa.gz

(time $bwa mem hg19_GL.fa.gz ../fakeq/hg19_GL.1.fastq.gz | samtools view -Sub - |
$samtools sort - hg19_GL.1.bwa) &> bwa_mem.log && 
$samtools index hg19_GL.1.bwa.bam
$samtools view -H hg19_GL.1.bwa.bam > hg19_GL.1.hdr.sam

alias revert="mv hg19_GL.1.uqct_prev.bam hg19_GL.1.uqct.bam && mv hg19_GL.1.uqct_prev.bam.bai hg19_GL.1.uqct.bam.bai"

../tst hg19_GL.fa.gz -c -m 12 &&
mv hg19_GL.1.uqct.bam hg19_GL.1.uqct_prev.bam &&
mv hg19_GL.1.uqct.bam.bai hg19_GL.1.uqct_prev.bam.bai &&
(cat hg19_GL.1.hdr.sam
../uqct ../fakeq/hg19_GL.1.fastq.gz ../hg19_GL.2b.gz -l 51) |
samtools view -Sub - |
$samtools sort - hg19_GL.1.uqct &&
$samtools index hg19_GL.1.uqct.bam &&
~/dev/git/IGV/igv.sh -b batchfile


cat << EOF > batchfile
new
load hg19_GL.1.bwa.bam
load hg19_GL.1.uqct_prev.bam
load hg19_GL.1.uqct.bam
goto GL000207.1
EOF


diff -u <(samtools view hg19_GL.1.uqct.bam)\
    <(samtools view hg19_GL.1.uqct_prev.bam) | gview -

#########################
make clean && DEFINES="-DKEY_LENGTH=17" make

ln -s /net/NGSanalysis/ref/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa
./uqct Homo_sapiens.GRCh38.dna.primary_assembly.fa -l 51


(cat ../hg38_hdr.sam
../uqct ../realfq/3387_1_MH_K27_05_CAATGGAA_L003_R1_001.fastq.gz ../Homo_sapiens.GRCh38.dna.primary_assembly.2b.gz -l 51) |
samtools view -Sub - |
$samtools sort - 3387_1_hg38.uqct &&
$samtools index 3387_1_hg38.uqct.bam

$samtools view -H hg19_GL.1.bwa.bam > hg19_GL.1.hdr.sam



