#!/bin/bash

DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
GAP2SEQ=$DIR/Gap2Seq
GAPCUTTER=$DIR/GapCutter
GAPMERGER=$DIR/GapMerger

usage() {
  echo "Gap2Seq [parameters]" 1>&2
  echo "" 1>&2
  echo "Required parameters:" 1>&2
  echo "-scaffolds <FASTA/Q file>    scaffolds to be gap filled" 1>&2
  echo "-filled <FASTA file>         output file for filled scaffolds" 1>&2
  echo "-reads <FASTA/Q files>       short reads, several files can be specified as a list separated by ','" 1>&2
  echo "" 1>&2
  echo "Optional parameters:" 1>&2
  echo "-max-mem <float>             maximum memory usage of DP table computation in gigabytes (excluding DBG) [default 20]" 1>&2
  echo "-fuz <int>                   number of nucleotides to ignore on gap fringes  [default 10]" 1>&2
  echo "-dist-error <int>            maximum error in gap estimates  [default 500]" 1>&2
  echo "-solid <int>                 threshold for solid k-mers for building the DBG [default 2]" 1>&2
  echo "-k                           kmer length for DBG  [default 31]" 1>&2
  echo "-nb-cores                    number of cores to use [default 0 (all cores)]" 1>&2
  echo "-verbose                     verbosity level (currently does not affect much?)  [default 1]" 1>&2
  echo "-help                        display help about possible options" 1>&2
  echo "-all-upper                   If specified, all filled bases are in upper case." 1>&2
  echo "-unique                      If specified, only gaps with a unique path of best length are filled." 1>&2
  exit 1
}

scaffolds=''
filled=''
reads=''
k=31
fuz=10
solid=2
dist=500
nb=0
verbose=1
maxmem=20
allupper=0
unique=0

while true; do
  case $1 in
    -scaffolds)
      shift
      scaffolds=$1 ;;

    -filled)
      shift
      filled=$1 ;;

    -reads)
      shift
      reads=$1 ;;

    -k)
      shift
      k=$1 ;;

    -dist-error)
      shift
      dist=$1 ;;

    -fuz)
      shift
      fuz=$1 ;;

    -solid)
      shift
      solid=$1 ;;

    -max-mem)
      shift
      maxmem=$1 ;;

    -nb-cores)
      shift
      nb=$1 ;;

    -verbose)
      shift
      verbose=$1 ;;

    -all-upper)
      shift
      allupper=1 ;;

    -unique)
      shift
      unique=1 ;;

    -*)
      echo "$0: Unrecognized option $1" >&2
      usage ;;

    *)
      if [[ $1 = '' ]]; then
        break
      fi

      shift ;;
  esac
done

if [[ $scaffolds = '' ]]; then
  echo "Error : Option '-scaffolds' is mandatory" >&2
  usage
fi

if [[ $filled = '' ]]; then
  echo "Error : Option '-filled' is mandatory" >&2
  usage
fi

if [[ $reads = '' ]]; then
  echo "Error : Option '-reads' is mandatory" >&2
  usage
fi

gaps=tmp.gaps
contigs=tmp.contigs
gapfilled=tmp.fill

options="-scaffolds $gaps -filled $gapfilled -reads $reads -k $k -fuz $fuz -nb-cores $nb -solid $solid -max-mem $maxmem -dist-error $dist -verbose $verbose"

if [[ $allupper = 1 ]]; then
  options="$options -all-upper"
fi

if [[ $unique = 1 ]]; then
  options="$options -unique"
fi

$GAPCUTTER -scaffolds $scaffolds -gaps $gaps -contigs $contigs -k $k -fuz $fuz || exit 1
$GAP2SEQ $options || exit 1
$GAPMERGER -scaffolds $filled -contigs $contigs -gaps $gapfilled || exit 1

rm $gaps $contigs $gapfilled
