#!/bin/bash
source $(dirname $0)/../environment.sh

cd $hdesdir/runtests
LASTKW=
ls -1 *.fa | sed -n -r 's/^(KW([0-9]+)_RL([0-9]+)_no[0-9]+)\.fa$/\1 \2 \3/p' |
while read BN KW RL; do
  if [ "$KW" != "$LASTKW" ]; then
    cd $hdesdir;
    make clean
    DEFINES="-DKEY_LENGTH=${KW}" make;
    cd -
    LASTKW="$KW"
  fi
  rm $BN.{2b,nn,bd,ub,kc,uq};
  valgrind $hdesdir/uqct ${BN}.fa -l ${RL}
  valgrind $hdesdir/uqct ${BN}.fq ${BN}.2b -l ${RL} |
  $samtools view -Sub - |
  $samtools sort -m 20000000 - ${BN} &&
  $samtools index ${BN}.bam

done
