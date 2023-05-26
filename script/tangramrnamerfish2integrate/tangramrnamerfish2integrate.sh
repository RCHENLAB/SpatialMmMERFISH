#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

args=$(getopt -o hd:b:r:l:s:e:g:m:c -l help,outdir:,bname:,rna:,clusterlabel:,seed:,condaenv:,gpu:,mode:,cpu --name "$0" -- "$@") || exit "$?"
eval set -- "$args"

absdir=$(dirname $(readlink -f "$0"))
scriptname=$(basename "$0" .sh)

outdir=.
clusterlabel=celltype
seed=12345
mode=cells
device=cuda
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
		-r|--rna)
			rna=$2
			shift 2
			;;
		-l|--clusterlabel)
			clusterlabel=$2
			shift 2
			;;
		-s|--seed)
			seed=$2
			shift 2
			;;
		-e|--condaenv)
			condaenv=$2
			shift 2
			;;
		-g|--gpu)
			gpu=$2
			device=cuda
			shift 2
			;;
		-m|--mode)
			mode=$2
			shift 2
			;;
		-c|--cpu)
			device=cpu
			shift
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
[[ $rna ]] || { echo "$scriptname.sh: -r|--rna must be specified!"; exit 1; }

function cmd {
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

if ! [[ $gpu ]]
then
	host=$(hostname)
	if [[ $host == mhgcp-g00.grid.bcm.edu ]]
	then
		gpu=$((RANDOM%2))
	elif [[ $host == mhgcp-g01.grid.bcm.edu ]]
	then
		gpu=$((RANDOM%4))
	else
		gpu=$((RANDOM%$(ngup.sh)))
	fi
fi
export CUDA_VISIBLE_DEVICES=$gpu

local f=$1
set -xe
hdf5ls.sh "$f" "$rna"
mkdir -p "$outdir" && pycmd.sh \
	-e "rna='$rna'" \
	-e "merfish='$f'" \
	-e "outdir='$outdir'" \
	-e "bname='$bname'" \
	-e "clusterlabel='$clusterlabel'" \
	-e "seed=$seed" \
	-e "mode='$mode'" \
	-e "device='$device'" \
	-s "$absdir/python/$scriptname.py"
}

if (($#))
then
	cmd "$@"
fi
