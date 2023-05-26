#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

args=$(getopt -o hd:b:r:g:c:A:B:H:W:t:e: -l help,outdir:,bname:,rowmean:,group:,cond:,condA:,condB:,width:,height:,numthreads:,condaenv: --name "$0" -- "$@") || exit "$?"
eval set -- "$args"

absdir=$(dirname $(readlink -f "$0"))
scriptname=$(basename "$0" .sh)

outdir=.
rowmean=-1
group=subtype
cond=layer
condA=GCL
condB=INL
width=8
height=6
numthreads=4
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
		-r|--rowmean)
			rowmean=$2
			shift 2
			;;
		-g|--group)
			group=$2
			shift 2
			;;
		-c|--cond)
			cond=$2
			shift 2
			;;
		-A|--condA)
			condA=$2
			shift 2
			;;
		-B|--condB)
			condB=$2
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

source traptimes
SECONDS=0
local f=$1
set -ex
hdf5ls.sh "$f"
mkdir -p "$outdir" && Rvanilla.sh \
	-e "infile='$f'" \
	-e "outdir='$outdir'" \
	-e "bname='$bname'" \
	-e "rowmean=$rowmean" \
	-e "group='$group'" \
	-e "cond='$cond'" \
	-e "condA='$condA'" \
	-e "condB='$condB'" \
	-e "width=$width" \
	-e "height=$height" \
	-e "numthreads=$numthreads" \
	-e "source('$absdir/R/$scriptname.R')"
set +x
echo "Time elapsed: $(date -u -d @$SECONDS +'%H:%M:%S')"
}

if (($#))
then
	cmd "$@"
fi
