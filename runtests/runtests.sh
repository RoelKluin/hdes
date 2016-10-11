#!/bin/bash
source $(dirname -- "$(readlink -f $0)")/../environment.sh
echo $hdesdir
GDB_OPT="-O0" # -Og should be less slow but some still get optimized away.

cd $hdesdir

mode=valgrind
[ "$1" = "-g" ] && mode=gdb
[ "$1" = "-g2" ] && mode=gdb2

FASTAS=($(ls -1 runtests/*.fa))
for i in $(seq 0 $((${#FASTAS[@]}-1))); do
    echo -e "$(echo -n $i | perl -n -e 'print chr($_+48)'))\t${FASTAS[${i}]}"
done
echo
read -n 1 -p "select file [q]uit or all (return)" ask;
if [ "$ask" == "q" ]; then
   exit 0;
else
    ask=$(echo -n ${ask}| perl -n -e 'print (ord($_)-48)')
    if [[ $ask -le ${#FASTAS[@]} ]]; then
        FASTAS=("${FASTAS[$ask]}")
    fi
fi
echo

make cleartest
LASTKW=
echo "${FASTAS[@]}" | tr " " "\n" | sed -n -r 's/^(runtests\/KW([0-9]+)_RL([0-9]+)_no[0-9]+)\.fa$/\1 \2 \3/p' |
while read BASE KW RL; do
  echo "--------------------------[ $BASE $KW $RL ]--------------------------";
  if [ "$KW" != "$LASTKW" ]; then
    make clean
    [ "$mode" = "gdb" -o "$mode" = "gdb2" ] && DEFINES="-DKEY_LENGTH=${KW}" OPT="$GDB_OPT -ggdb3" make || DEFINES="-DKEY_LENGTH=${KW}" make;
    LASTKW="$KW"
  fi
  if [ "$mode" = "gdb" ];then
    echo "DEFINES='-DKEY_LENGTH=${KW}' gdb --command=gdb_index --args ./uqct -f ${BASE}.fa -l ${RL}"
    xterm -geometry 119x99-0+0 -b 0 -e \
        "DEFINES='-DKEY_LENGTH=${KW}' gdb --command=gdb_index --args ./uqct -f ${BASE}.fa -l ${RL}"
#    xterm -e "gdb --args uqct -f ${BASE}.fq ${BASE}.2b -l ${RL}"
  elif [ "$mode" = "gdb2" ];then
    valgrind --leak-check=full ./uqct ${BASE}.fa -l ${RL}
    echo "DEFINES='-DKEY_LENGTH=${KW}' gdb --command=gdb_map --args ./uqct -f ${BASE}.fq ${BASE}.2b -l ${RL}"
    xterm -geometry 119x99-0+0 -b 0 -e \
        "DEFINES='-DKEY_LENGTH=${KW}' gdb --command=gdb_map --args ./uqct ${BASE}.fq ${BASE}.2b -l ${RL}"
  else
    valgrind --leak-check=full ./uqct ${BASE}.fa -l ${RL}
    valgrind --leak-check=full ./uqct ${BASE}.fq ${BASE}.2b -l ${RL} |
    $samtools view -Sub - |
    $samtools sort -m 20000000 - ${BASE} &&
    $samtools index ${BASE}.bam
  fi
done
echo
echo
