#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

args=$(getopt -o hd:b:H:W:ne: -l help,outdir:,bname:,height:,width:,norm,condaenv: --name "$0" -- "$@") || exit "$?"
eval set -- "$args"

absdir=$(dirname $(readlink -f "$0"))
scriptname=$(basename "$0" .sh)

outdir=.
height=5
width=5
norm=False
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
		-H|--height)
			height=$2
			shift 2
			;;
		-W|--width)
			width=$2
			shift 2
			;;
		-n|--norm)
			norm=True
			shift
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

local f=$1
set -x
mkdir -p "$outdir" && pycmd.sh \
	-e "f='$f'" \
	-e "outdir='$outdir'" \
	-e "bname='$bname'" \
	-e "width=$width" \
	-e "height=$height" \
	-e "norm=$norm" \
	-s "$absdir/python/$scriptname.py"
}

if (($#))
then
	cmd "$@"
fi
