#!/bin/bash
source $(dirname -- "$(readlink -f $0)")/../environment.sh
echo $hdesdir

cd $hdesdir/runtests

FASTAS=($(ls -1 *.fa))
for i in $(seq 0 $((${#FASTAS[@]}-1))); do
  echo -e "$i)\t${FASTAS[${i}]}"
done
echo
read -n 1 -p "select file [q]uit or all (return)" ask;
if [ "$ask" == "q" ]; then
   exit 0;
elif [[ "$ask" =~ ^[0-9]+$ ]]; then
  FASTAS=("${FASTAS[${ask}]}")
fi
echo

LASTKW=
echo "${FASTAS[@]}" | tr " " "\n" | sed -n -r 's/^(KW([0-9]+)_RL([0-9]+)_no[0-9]+)\.fa$/\1 \2 \3/p' |
while read BN KW RL; do
  echo "--------------------------[ $BN $KW $RL ]--------------------------";
  if [ "$KW" != "$LASTKW" ]; then
    cd $hdesdir;
    make clean
    DEFINES="-DKEY_LENGTH=${KW}" make;
    cd -
    LASTKW="$KW"
  fi
  rm $BN.{2b,nn,bd,ub,kc,uq};
  valgrind --leak-check=full $hdesdir/uqct ${BN}.fa -l ${RL}
  valgrind --leak-check=full $hdesdir/uqct ${BN}.fq ${BN}.2b -l ${RL} |
  $samtools view -Sub - |
  $samtools sort -m 20000000 - ${BN} &&
  $samtools index ${BN}.bam

done
echo
echo
