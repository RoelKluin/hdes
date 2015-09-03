#!/bin/bash
cd /home/roel/dev/git/hdes/
die() {
  echo -e "$1"
  [ -z "$2" ] && exit 1 || exit $2;
}
USAGE="$0 [-c|--clean] [-p|--part] [-m|--make <KEY_LENGTH>] [-L|--length <READLENGTH>] [-v|--valgrind] [-l|--leak-check] fasta.gz"

TEMP=`getopt -o m:L:pcvl --long make:,length:,clean,part,valgrind,leak-check -n "$0" -- "$@"`
if [ $? != 0 ]; then echo "Terminating..." >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

PART=
MAKE=
CLEAN=
VALGRIND=
LENGTH=51
while true; do
  case "$1" in
    -m|--make) MAKE="$2"; shift;;
    -L|--length) LENGTH="$2"; shift;;
    -p|--part) PART=1;;
    -c|--clean) CLEAN=1;;
    -v|--valgrind) [ -z "$VALGRIND" ] && VALGRIND=valgrind;;
    -l|--leak-check) VALGRIND="valgrind --leak-check=full";;
    -- ) shift; break;;
    * ) break;;
  esac
  shift;
done
[ $# -ne 1 ] && die "$USAGE\n\nneed fasta.gz"
F="$1"
[[ "$F" =~ \.gz$ ]] || die "$USAGE\n\nnot gzipped"
[ -e "$F" ] || "$USAGE\n\nno such file: $F"


BN="$(basename "$F" ".gz")"
[ -e "${BN}.err" ] && mv "${BN}.err" "${BN}.old.err"

if [ -n "$PART" ]; then
  for f in "${BN%.fa}".{ub,uq}.gz; do
    mv $f ${f%.gz}.old.gz 2> /dev/null;
  done
elif [ -n "$CLEAN" ]; then
  for f in "${BN%.fa}".{2b,nn,bd,ub,kc,uq}.gz; do
    mv $f ${f%.gz}.old.gz 2> /dev/null;
  done
  make clean || exit 1;
fi
if [ -n "$MAKE" ]; then
  DEFINES="-DKEY_LENGTH=$MAKE" make 2>&1 | tee ${BN%.fa}.err || exit 1
fi
$VALGRIND ./uqct "$F" -l $LENGTH 2>&1 | tee -a ${BN%.fa}.err || exit 1

