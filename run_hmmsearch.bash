#!/bin/bash
set -euo pipefail
usage() { echo "Usage: $0 [-f families] [-t threads] [-s suffix] input_fasta" 1>&2; exit 1; }
tmpdir=""
cleanup() { if [[ ! -z "$tmpdir" && -e $tmpdir ]]; then echo "cleaning up"; rm -r $tmpdir; fi }
trap cleanup EXIT

families=legfed_v1_0
threads=64
suffix=""
while getopts ":t:s:f:" o; do
    case "${o}" in
        f)
            families=${OPTARG}
            ;;
        t)
            threads=${OPTARG}
            ;;
        s)
            suffix=${OPTARG}
            ;;
        h|*)
            usage
            ;;
    esac
done
shift $((OPTIND-1))
if (( $# != 1 )); then usage; fi
input=$1
hmm=""
case $families in
	legume.fam3)
            hmm=$DATA/legume.fam3.hmm
	    ;;
	legfed_v1_0)
            hmm=$DATA/legfed_v1_0.hmm
	    ;;
	phytozome_10_2)
            hmm=$DATA/Angiosperm.hmm
	    ;;
	PANTHER_17_0)
            hmm=$DATA/PANTHER17.0_hmm
	    ;;
esac
if [[ -z $hmm ]]; then
    echo "no hmm known for families = $families"
    exit 1
fi
#hmmsearch doesn't handle compressed fasta natively
if [[ $1 =~ \.gz$ ]]; then 
    tmpprefix=`basename $0`
    tmpdir=`mktemp -d ${tmpprefix}.XXXXXX`
    tmpname=`basename $1`
    tmp=${tmpdir}/${tmpname/.gz/}
    echo "uncompressing input to $tmp"
    gzip -dc $1 > $tmp
    input=$tmp
fi
dir=hmmsearch_${families}${suffix}
mkdir -p $dir
hmmsearch --cpu $threads $hmm $input | gzip -c > ${dir}/proteins.hmmsearch.tbl.gz
