#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

args=$(getopt -o hd:b:x:y:e:g:l:s:H:W:n: -l help,outdir:,bname:,xcoord:,ycoord:,layer:,group:,label:,split:,width:,height:,nreplicate: --name "$0" -- "$@") || exit "$?"
eval set -- "$args"

absdir=$(dirname $(readlink -f "$0"))
scriptname=$(basename "$0" .sh)

outdir=.
xcoord=final_x
ycoord=final_y
layer=layer
group=GCL
label=subtype
split=sectionid
width=6
height=5
nreplicate=1000
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
		-x|--xcoord)
			xcoord=$2
			shift 2
			;;
		-y|--ycoord)
			ycoord=$2
			shift 2
			;;
		-e|--layer)
			layer=$2
			shift 2
			;;
		-g|--group)
			group=$2
			shift 2
			;;
		-l|--label)
			label=$2
			shift 2
			;;
		-s|--split)
			split=$2
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
		-n|--nreplicate)
			nreplicate=$2
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
source traptimes
SECONDS=0
local f=$1
set -ex
mkdir -p "$outdir" && Rvanilla.sh \
	-e "f='$f'" \
	-e "outdir='$outdir'" \
	-e "bname='$bname'" \
	-e "xcoord='$xcoord'" \
	-e "ycoord='$ycoord'" \
	-e "layer='$layer'" \
	-e "group='$group'" \
	-e "label='$label'" \
	-e "split='$split'" \
	-e "width=$width" \
	-e "height=$height" \
	-e "nreplicate=$nreplicate" \
	-e "source('$absdir/R/$scriptname.R')"
set +x
echo "Time elapsed: $(date -u -d @$SECONDS +'%H:%M:%S')"
}

if (($#))
then
	cmd "$@"
fi
