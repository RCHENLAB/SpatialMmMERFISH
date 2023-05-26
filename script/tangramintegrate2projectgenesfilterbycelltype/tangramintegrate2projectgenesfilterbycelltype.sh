#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

args=$(getopt -o hd:b:m:r:n:t:w:e: -l help,outdir:,bname:,map:,rna:,rnatype:,merfishtype:,totalweight:,condaenv: --name "$0" -- "$@") || exit "$?"
eval set -- "$args"

absdir=$(dirname $(readlink -f "$0"))
scriptname=$(basename "$0" .sh)

outdir=.
totalweight=-1
while true
do
	case "$1" in
		-h|--help)
			exec cat "$absdir/$scriptname.txt"
			;;
		-d|--outdir)
			outdir=$2
			shift 2
			;;
		-b|--bname)
			bname=$2
			shift 2
			;;
		-m|--map)
			map=$2
			shift 2
			;;
		-r|--rna)
			rna=$2
			shift 2
			;;
		-n|--rnatype)
			rnatype=$2
			shift 2
			;;
		-t|--merfishtype)
			merfishtype=$2
			shift 2
			;;
		-w|--totalweight)
			totalweight=$2
			shift 2
			;;
		-e|--condaenv)
			condaenv=$2
			shift 2
			;;
		--)
			shift
			break
			;;
		*)
			echo "$0: not implemented option: $1" >&2
			exit 1
			;;
	esac
done

[[ $bname ]] || { echo "$scriptname.sh: -b|--bname must be specified!"; exit 1; }
[[ $map ]] || { echo "$scriptname.sh: -m|--map must be specified!"; exit 1; }
[[ $rna ]] || { echo "$scriptname.sh: -r|--rna must be specified!"; exit 1; }
[[ $rnatype ]] || { echo "$scriptname.sh: -n|--rnatype must be specified!"; exit 1; }
[[ $merfishtype ]] || { echo "$scriptname.sh: -t|--merfishtype must be specified!"; exit 1; }

if [[ $condaenv ]]
then
	if [[ $MAMBA_ROOT_PREFIX ]]
	then
		source "$MAMBA_ROOT_PREFIX/etc/profile.d/micromamba.sh"
		micromamba activate "$condaenv"
	else
		source "$(conda info --base)/etc/profile.d/conda.sh"
		conda activate "$condaenv"
	fi
fi

set -xe
hdf5ls.sh "$map" "$rna"
mkdir -p "$outdir" && pycmd.sh \
	-e "map='$map'" \
	-e "rna='$rna'" \
	-e "outdir='$outdir'" \
	-e "bname='$bname'" \
	-e "rnatype='$rnatype'" \
	-e "merfishtype='$merfishtype'" \
	-e "totalweight=$totalweight" \
	-s "$absdir/python/$scriptname.py"
