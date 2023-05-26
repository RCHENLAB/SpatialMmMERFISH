#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import sys
from exeR.exeR import exeR
from pathlib import Path
from multiprocess import Pool
import click

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(resolve_path=True), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', multiple=True, type=click.STRING, help='Bnames for infile. Default: basenames of infiles.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-H', '--height', type=click.FLOAT, default=5, show_default=True, help='Height.')
@click.option('-W', '--width', type=click.FLOAT, default=6, show_default=True, help='Width.')
@click.option('-t', '--numthreads', type=click.INT, default=4, show_default=True, help='Number of threads.')
@click.option('-n', '--ntop', type=click.INT, default=15, show_default=True, help='Number of top GO terms.')
@click.option('-f', '--fontsize', type=click.FLOAT, default=10, show_default=True, help='Font size.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True), nargs=-1)
def main(outdir, bname, condaenv, height, width, numthreads, ntop, fontsize, infile):
	"""
Plot the barplot and dotplot of enriched GO terms.

INFILE is a enrichGO() RDS file[s].

\b
Example:
  indir=$(mrrdir.sh ../deseq2enrichgo)
  outdir=$(mrrdir.sh)
  slurmtaco.sh -t 2 -m 20G -- scrnaenrichgo2dotplot -d "$outdir" -- "$indir"/*.rds

\b
See also:
  Upstream:
    deseq2enrichgo.sh
    enrichgobygeneset.sh
  Related:
    scrnaenrichgo2enrichplot.sh
  Depends:
    R/clusterProfiler
    R/enrichplot

\b
Date: 2023/04/11
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	if len(bname)!=len(infile):
		bname=[Path(file).stem for file in infile]

	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f'{absdir}/R/{scriptname}.R'
	Path(outdir).mkdir(parents=True, exist_ok=True)
	# os.chdir(outdir)

	def cmd(*args):
		name, file=args
		exprs=[
			f"outdir='{outdir}'",
			f"bname='{name}'",
			f"height={height}",
			f"width={width}",
			f"ntop={ntop}",
			f"fontsize={fontsize}",
			f"infile='{file}'",
			]
		exeR.callback(exprs, script=script, condaenv=condaenv, verbose=True)

	with Pool(processes=numthreads) as p:
		p.starmap(cmd, zip(bname, infile))

	return 0

if __name__ == "__main__":
	sys.exit(main())
