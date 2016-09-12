#!/bin/bash

source $(dirname $0)/environment.sh


cd $hdesdir
die() {
  echo -e "$1"
  [ -z "$2" ] && exit 1 || exit $2;
}
USAGE="$0 [-g] [-c|--clean] [-C|--callgrind] [-p|--part] [-m|--make <KEY_LENGTH>] [-L|--length <READLENGTH>] [-v|--valgrind] [-l|--leak-check] fasta.gz"

TEMP=`getopt -o m:L:pcvlg --long make:,length:,clean,callgrind,part,valgrind,leak-check:gdb -n "$0" -- "$@"`
if [ $? != 0 ]; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

PART=
MAKE=
CLEAN=
PRE=
LENGTH=51
while true; do
  case "$1" in
    -m|--make) MAKE="$2"; shift;;
    -L|--length) LENGTH="$2"; shift;;
    -p|--part) PART=1;;
    -g|--gdb) PRE="gdb --args";;
    -c|--clean) CLEAN=1;;
    -C|--callgrind) [ -z "$PRE" ] && PRE="valgrind --tool=callgrind";;
    -v|--valgrind) [ -z "$PRE" ] && PRE=valgrind;;
    -l|--leak-check) PRE="valgrind --leak-check=full";;
    -- ) shift; break;;
    * ) break;;
  esac
  shift;
done
[ $# -ne 1 ] && die "$USAGE\n\nneed fasta.gz"
F="$1"
[ -e "$F" ] || "$USAGE\n\nno such file: $F"

BN="${F%.fa*}"



[ -e "${BN}.err" ] && mv "${BN}.err" "${BN}.old.err"

if [ -n "$PART" ]; then
  for x in ub uq; do
    mv ${F/.fa/.${x}} ${F/.fa/.${x}.old};
  done
elif [ -n "$CLEAN" ]; then
  for x in 2b nn bd ub kc uq; do
    mv ${F/.fa/.${x}} ${F/.fa/.${x}.old};
  done
  make clean || exit 1;
fi

if [ -n "$MAKE" ]; then
  DEFINES="-DKEY_LENGTH=$MAKE" make 2>&1 | tee ${BN}.err || exit 1
fi
echo -e "Running\n$PRE ./uqct \"$F\" -l $LENGTH 2>&1 | tee -a ${BN}.err"
$PRE ./uqct "$F" -l $LENGTH 2>&1 | tee -a ${BN}.err || exit 1

