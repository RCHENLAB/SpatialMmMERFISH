#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

args=$(getopt -o hd:b:l:W:H:n: -l help,outdir:,bname:,label:,width:,height:,ntop: --name "$0" -- "$@") || exit "$?"
eval set -- "$args"

absdir=$(dirname $(readlink -f "$0"))
scriptname=$(basename "$0" .sh)

outdir=.
width=6
height=5
ntop=2000
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
		-l|--label)
			label=$2
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

[[ $bname ]] || { echo "$scriptname.sh: -b|--bname must be specified!"; exit 1; }
[[ $label ]] || { echo "$scriptname.sh: -l|--label must be specified!"; exit 1; }

function cmd {
source traptimes
SECONDS=0
local f=$1
set -x
mkdir -p "$outdir" && Rvanilla.sh \
	-e "f='$f'" \
	-e "outdir='$outdir'" \
	-e "bname='$bname'" \
	-e "label='$label'" \
	-e "width=$width" \
	-e "height=$height" \
	-e "ntop=$ntop" \
	-e "source('$absdir/R/$scriptname.R')"
set +x
echo "Time elapsed: $(date -u -d @$SECONDS +'%H:%M:%S')"
}

if (($#))
then
	cmd "$@"
fi
