#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

args=$(getopt -o ho:W:H:b:x:m:T -l help,outfile:,width:,height:,breaks:,xlab:,main:,transform,density,fixedwidth,lower:,higher:,leftclose,combine --name "$0" -- "$@") || exit "$?"
eval set -- "$args"

absdir=$(dirname $(readlink -f "${BASH_SOURCE[0]}"))
scriptname=$(basename "${BASH_SOURCE[0]}" .sh)

width=5
height=5
breaks=20
xlab=x
main=
transform=F
density=F
fixedwidth=F
right=T
combine=F
while true
do
	case "$1" in
		-h|--help)
			exec cat "$absdir/$scriptname.txt"
			;;
		-o|--outfile)
			outfile=$2
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
		-b|--breaks)
			breaks=$2
			shift 2
			;;
		-x|--xlab)
			xlab=$2
			shift 2
			;;
		-m|--main)
			main=$2
			shift 2
			;;
		-T|--transform)
			transform=T
			shift
			;;
		--density)
			density=T
			shift
			;;
		--fixedwidth)
			fixedwidth=T
			shift
			;;
		--lower)
			lower=$2
			shift 2
			;;
		--higher)
			higher=$2
			shift 2
			;;
		--leftclose)
			right=F
			shift
			;;
		--combine)
			combine=T
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

[ "$outfile" ] || { echo "$scriptname.sh: -o|--outfile must be specified!"; exit 1; }

if [[ $fixedwidth == T ]]
then
	[ "$lower" ] || { echo "$scriptname.sh: --lower must be specified for fixed width!"; exit 1; }
	[ "$higher" ] || { echo "$scriptname.sh: --higher must be specified for fixed width!"; exit 1; }
fi

mkdir -p "$(dirname "$outfile")" && exec R --slave --no-save --no-restore --no-init-file \
	-e "width=$width" \
	-e "height=$height" \
	-e "outfile='$outfile'" \
	-e "breaks=$breaks" \
	-e "xlab='$xlab'" \
	-e "main='$main'" \
	-e "transform='$transform'" \
	-e "density=$density" \
	-e "fixedwidth=$fixedwidth" \
	-e "lower=$lower" \
	-e "higher=$higher" \
	-e "right=$right" \
	-e "combine=$combine" \
	-e "source('$absdir/R/$scriptname.R')"
