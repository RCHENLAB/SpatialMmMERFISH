#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

args=$(getopt -o hd:b:g:e:m:H:W:c:snf: -l help,outdir:,bname:,genefile:,condaenv:,method:,height:,width:,colorkey:,shownumber,normonly,scatter,filter: --name "$0" -- "$@") || exit "$?"
eval set -- "$args"

absdir=$(dirname $(readlink -f "$0"))
scriptname=$(basename "$0" .sh)

outdir=.
method=pearson
width=5
height=4
colorkey=0.5
shownumber=F
normonly=F
scatter=F
vfilter=NULL
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
		-g|--genefile)
			genefile=$2
			shift 2
			;;
		-e|--condaenv)
			condaenv=$2
			shift 2
			;;
		-m|--method)
			method=$2
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
		-c|--colorkey)
			colorkey=$2
			shift 2
			;;
		-s|--shownumber)
			shownumber=T
			shift
			;;
		-n|--normonly)
			normonly=T
			shift
			;;
		--scatter)
			scatter=T
			shift
			;;
		-f|--filter)
			vfilter=$2
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
[[ $genefile ]] || { echo "$scriptname.sh: -g|--genefile must be specified!"; exit 1; }

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
set -xe

hdf5ls.sh "$f"
mkdir -p "$outdir" && exec Rvanilla.sh \
	-e "infile='$f'" \
	-e "outdir='$outdir'" \
	-e "bname='$bname'" \
	-e "genefile='$genefile'" \
	-e "method='$method'" \
	-e "width=$width" \
	-e "height=$height" \
	-e "colorkey=$colorkey" \
	-e "shownumber=$shownumber" \
	-e "normonly=$normonly" \
	-e "scatter=$scatter" \
	-e "vfilter=$vfilter" \
	-e "source('$absdir/R/$scriptname.R')"

set +x
echo "Time elapsed: $(date -u -d @$SECONDS +'%H:%M:%S')"
}

if (($#))
then
	cmd "$@"
fi
