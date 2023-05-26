#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

args=$(getopt -o hd:q:l:t:s: -l help,outdir:,qthr:,logfcthr:,numthreads:,species: --name "$0" -- "$@") || exit "$?"
eval set -- "$args"

absdir=$(dirname $(readlink -f "$0"))
scriptname=$(basename "$0" .sh)

outdir=.
numthreads=$(nproc)
qthr=0.05
logfcthr=1
species=hs
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
		-q|--qthr)
			qthr=$2
			shift 2
			;;
		-l|--logfcthr)
			logfcthr=$2
			shift 2
			;;
		-t|--numthreads)
			numthreads=$2
			shift 2
			;;
		-s|--species)
			species=$2
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
local bname=$(basename.sh -s .txt -s .txt.gz -s .tsv -s .tsv.gz -s .rnk -s .rnk.gz -- "$f")
local symbolup=$outdir/${bname}_symbol_up.txt.gz
local symboldown=$outdir/${bname}_symbol_down.txt.gz
local outup=$outdir/${bname}_GO_up.txt.gz
local outdown=$outdir/${bname}_GO_down.txt.gz
set -x
cat.sh -q "$f" | hcut.sh padj symbol log2FoldChange | tail -n +2 | tsvfilterbycol.sh --le "$qthr" | tsvfilterbycol.sh -c 3 --ge "$logfcthr" | cut -f 2 | tofile.sh -o "$symbolup"
cat.sh -q "$f" | hcut.sh padj symbol log2FoldChange | tail -n +2 | tsvfilterbycol.sh --le "$qthr" | tsvfilterbycol.sh -c 3 --le -- -"$logfcthr" | cut -f 2 | tofile.sh -o "$symboldown"

enrichgobygeneset.sh --"$species" -o "$outup" -- "$symbolup"
enrichgobygeneset.sh --"$species" -o "$outdown" -- "$symboldown"
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
