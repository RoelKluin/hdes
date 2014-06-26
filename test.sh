
#make && zgrep -E -o "^[ACTGN]{6}" /home/roel/dev/rproject/rmap/hg19.fa.gz | head -n 800000 | valgrind ./uqct 2>&1 | less

#time zgrep -E -o "^[ACTGN]{6}" /home/roel/dev/rproject/rmap/hg19.fa.gz | head -n 800000 | ./uqct


make && zcat /home/roel/dev/git/bwa/my/bwa/test/ce_med.fq.gz | valgrind ./uqct 1>/dev/null 2> valg.err

make && zcat /home/roel/dev/git/bwa/my/bwa/test/ce_med.fq.gz | ./uqct | less
