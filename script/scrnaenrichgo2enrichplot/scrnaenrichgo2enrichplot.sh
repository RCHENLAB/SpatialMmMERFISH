#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

args=$(getopt -o hd:W:H:t:n: -l help,outdir:,width:,height:,numthreads:,ntop: --name "$0" -- "$@") || exit "$?"
eval set -- "$args"

absdir=$(dirname $(readlink -f "$0"))
scriptname=$(basename "$0" .sh)

outdir=.
width=10
height=6
numthreads=$(nproc)
ntop=30
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
		-W|--width)
			width=$2
			shift 2
			;;
		-H|--height)
			height=$2
			shift 2
			;;
		-t|--numthreads)
			numthreads=$2
			shift 2
			;;
		-n|--ntop)
			ntop=$2
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

function cmd {
local f=$1
local bname=$(basename.sh -s .rds -- "$f")
set -x
mkdir -p "$outdir" && Rvanilla.sh \
	-e "infile='$f'" \
	-e "outdir='$outdir'" \
	-e "bname='$bname'" \
	-e "width=$width" \
	-e "height=$height" \
	-e "ntop=$ntop" \
	-e "source('$absdir/R/$scriptname.R')"
}

if (($#))
then
	if (($#==1))
	then
		cmd "$@"
	else
		source env_parallel.bash
		env_parallel -j "$numthreads" cmd ::: "$@"
	fi
fi
