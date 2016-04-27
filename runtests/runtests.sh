#!/bin/bash
source $(dirname -- "$(readlink -f $0)")/../environment.sh
echo $hdesdir
GDB_OPT="-O0" # -Og should be less slow but some still get optimized away.

cd $hdesdir

mode=valgrind
[ "$1" = "-g" ] && mode=gdb

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

make cleartest
LASTKW=
echo "${FASTAS[@]}" | tr " " "\n" | sed -n -r 's/^(KW([0-9]+)_RL([0-9]+)_no[0-9]+)\.fa$/\1 \2 \3/p' |
while read BN KW RL; do
  echo "--------------------------[ $BN $KW $RL ]--------------------------";
  if [ "$KW" != "$LASTKW" ]; then
    make clean
    [ "$mode" = "gdb" ] && DEFINES="-DKEY_LENGTH=${KW}" OPT="$GDB_OPT -ggdb3" make || DEFINES="-DKEY_LENGTH=${KW}" make;
    LASTKW="$KW"
  fi
  rm $BN.{2b,nn,bd,ub,kc,uq};
  if [ "$mode" = "gdb" ];then
    xterm -e "DEFINES='-DKEY_LENGTH=${KW}' gdb --args ./uqct -f runtests/${BN}.fa -l ${RL}"
#    xterm -e "gdb --args uqct -f runtests/${BN}.fq runtests/${BN}.2b -l ${RL}"
  else
    valgrind --leak-check=full ./uqct runtests/${BN}.fa -l ${RL}
    valgrind --leak-check=full ./uqct runtests/${BN}.fq runtests/${BN}.2b -l ${RL} |
    $samtools view -Sub - |
    $samtools sort -m 20000000 - runtests/${BN} &&
    $samtools index runtests/${BN}.bam
  fi
done
echo
echo
