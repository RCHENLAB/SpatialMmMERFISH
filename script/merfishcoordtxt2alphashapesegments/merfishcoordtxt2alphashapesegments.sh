#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

args=$(getopt -o hd:b:a:x:y:H:W: -l help,outdir:,bname:,alpha:,coordx:,coordy:,width:,height: --name "$0" -- "$@") || exit "$?"
eval set -- "$args"

absdir=$(dirname $(readlink -f "$0"))
scriptname=$(basename "$0" .sh)

outdir=.
alpha=50
coordx=center_x
coordy=center_y
width=8
height=6
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
		-a|--alpha)
			alpha=$2
			shift 2
			;;
		-x|--coordx)
			coordx=$2
			shift 2
			;;
		-y|--coordy)
			coordy=$2
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
local f=$1
set -x
mkdir -p "$outdir" && exec Rvanilla.sh \
	-e "f='$f'" \
	-e "outdir='$outdir'" \
	-e "bname='$bname'" \
	-e "alpha=$alpha" \
	-e "coordx='$coordx'" \
	-e "coordy='$coordy'" \
	-e "width=$width" \
	-e "height=$height" \
	-e "source('$absdir/R/$scriptname.R')"
}

if (($#))
then
	cmd "$@"
fi
