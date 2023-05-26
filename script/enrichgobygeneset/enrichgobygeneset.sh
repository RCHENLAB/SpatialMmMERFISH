#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

args=$(getopt -o ho:s:q: -l help,outfile:,species:,mm,hs,qthr: --name "$0" -- "$@") || exit "$?"
eval set -- "$args"

absdir=$(dirname $(readlink -f "$0"))
scriptname=$(basename "$0" .sh)

species=Hs
qthr=0.05
while true
do
	case $1 in
		-h|--help)
			exec cat "$absdir/${scriptname}.txt"
			;;
		-o|--outfile)
			outfile=$2
			shift 2
			;;
		-s|--species)
			species=$2
			shift 2
			;;
		--mm)
			species=Mm
			shift
			;;
		--hs)
			species=Hs
			shift
			;;
		-q|--qthr)
			qthr=$2
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

[[ $outfile ]] || { echo "$scriptname.sh: -o|--outfile must be specified!"; exit 1; }

function cmd {
local f=$1
local outdir=$(dirname "$outfile")
local bname=$(basename.sh -s .rds -s .txt -s .txt.gz -s .tsv -s .tsv.gz -- "$outfile")
set -x
mkdir -p "$outdir" && Rvanilla.sh \
	-e "infile='$f'" \
	-e "outdir='$outdir'" \
	-e "bname='$bname'" \
	-e "species='org.$species.eg.db'" \
	-e "qthr=$qthr" \
	-e "source('$absdir/R/$scriptname.R')"
}

if (($#))
then
	cmd "$@"
fi
