#!/bin/bash
if [ "$(hostname)" = "Z" ]; then
  samtools=/home/roel/dev/git/samtools/lh3/samtools/samtools
  bwa=/home/roel/dev/git/bwa/orig/bwa/bwa
  hdesdir=/home/roel/dev/git/hdes
else
  samtools=/net/NGSanalysis/apps/samtools/samtools-0.1.19/samtools
  bwa=/net/NGSanalysis/apps/bwa/bwa-0.7.12/bwa
  hdesdir=/net/NGSanalysis/dvl/roel/git/hdes
fi
