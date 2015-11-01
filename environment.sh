#!/bin/bash
if [ "$(hostname)" = "Z" ]; then
  samtools=/home/roel/dev/git/samtools/lh3/samtools/samtools
  bwa=/home/roel/dev/git/bwa/orig/bwa/bwa
  hdesdir=/home/roel/dev/git/hdes
  igv=/home/roel/dev/git/IGV/igv.sh
  java=/opt/java/jre1.8.0_60/bin/java
  artfq=/home/roel/dev/src/ArtificialFastqGenerator/ArtificialFastqGenerator.jar
else
  samtools=/net/NGSanalysis/apps/samtools/samtools-0.1.19/samtools
  bwa=/net/NGSanalysis/apps/bwa/bwa-0.7.12/bwa
  if [[  $(dirname -- "$(readlink -f $0)") =~ home_tmp ]]; then
    hdesdir=/net/NGSanalysis/dvl/roel/git/hdes/home_tmp
  else
    hdesdir=/net/NGSanalysis/dvl/roel/git/hdes
  fi
  igv=/home/roel/Downloads/igv/igv.sh
  java=/opt/java/jre1.8.0_45/bin/java
  artfq=/net/NGSanalysis/dvl/roel/git/hdes/external/ArtificialFastqGenerator.jar
fi
