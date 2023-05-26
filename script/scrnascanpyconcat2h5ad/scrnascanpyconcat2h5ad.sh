#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

args=$(getopt -o hd:b:j:t:rze: -l help,outdir:,bname:,join:,numthreads:,raw,zerofill,condaenv: --name "$0" -- "$@") || exit "$?"
eval set -- "$args"

absdir=$(dirname $(readlink -f "$0"))
scriptname=$(basename "$0" .sh)

outdir=.
join=inner
numthreads=4
raw=False
zerofill=False
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
		-j|--join)
			join=$2
			shift 2
			;;
		-t|--numthreads)
			numthreads=$2
			shift 2
			;;
		-r|--raw)
			raw=True
			shift
			;;
		-z|--zerofill)
			zerofill=True
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
	source "$(conda info --base)/etc/profile.d/conda.sh"
	conda activate "$condaenv"
fi

local f=$1
set -x
mkdir -p "$outdir" && exec "$absdir/python/$scriptname.py" "$f" "$outdir" "$bname" "$join" "$numthreads" "$raw" "$zerofill"
}

if (($#))
then
	cmd "$@"
fi
