#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

args=$(getopt -o hd:b:s:x:y:c: -l help,outdir:,bname:,shiftby:,coordx:,coordy:,stepx:,stepy:,finalx:,finaly:,ncol: --name "$0" -- "$@") || exit "$?"
eval set -- "$args"

absdir=$(dirname $(readlink -f "$0"))
scriptname=$(basename "$0" .sh)

outdir=.
shiftby=sectionid
coordx=center_x
coordy=center_y
stepx=4000
stepy=4000
finalx=final_x
finaly=final_y
ncol=5
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
		-s|--shiftby)
			shiftby=$2
			shift 2
			;;
		--coordx)
			coordx=$2
			shift 2
			;;
		--coordy)
			coordy=$2
			shift 2
			;;
		-x|--stepx)
			stepx=$2
			shift 2
			;;
		-y|--stepy)
			stepy=$2
			shift 2
			;;
		--finalx)
			finalx=$2
			shift 2
			;;
		--finaly)
			finaly=$2
			shift 2
			;;
		-c|--ncol)
			ncol=$2
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
set -xe
hdf5ls.sh "$f"
mkdir -p "$outdir" && pycmd.sh \
	-e "f='$f'" \
	-e "outdir='$outdir'" \
	-e "bname='$bname'" \
	-e "shiftby='$shiftby'" \
	-e "coordx='$coordx'" \
	-e "coordy='$coordy'" \
	-e "stepx=$stepx" \
	-e "stepy=$stepy" \
	-e "finalx='$finalx'" \
	-e "finaly='$finaly'" \
	-e "ncol=$ncol" \
	-s "$absdir/python/$scriptname.py"
}

if (($#))
then
	cmd "$@"
fi
