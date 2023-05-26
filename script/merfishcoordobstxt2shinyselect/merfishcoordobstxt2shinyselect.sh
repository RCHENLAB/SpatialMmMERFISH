#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

args=$(getopt -o hd:b: -l help,outdir:,bname: --name "$0" -- "$@") || exit "$?"
eval set -- "$args"

absdir=$(dirname $(readlink -f "$0"))
scriptname=$(basename "$0" .sh)

outdir=.
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
local f=$(abspath.sh "$1")
set -xe
rsync.sh -d "$outdir" -- "$absdir/template/main.Rmd"
m4 -D INFILE="$f" -D OUTFILE="$outdir/$bname.txt.gz" < "$absdir/template/main.R" > "$outdir/main.R"
rstudio.sh "$outdir/main.Rmd"
}

if (($#))
then
	cmd "$@"
fi
