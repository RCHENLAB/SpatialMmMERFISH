#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import sys
from exeR.exeR import exeR
from pathlib import Path
from multiprocess import Pool
import pandas as pd
import click

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(resolve_path=True), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', type=click.STRING, required=True, help='Bname.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-t', '--numthreads', type=click.INT, default=4, show_default=True, help='Number of threads.')
@click.option('-s', '--species', type=click.Choice(['hs', 'mm']), is_flag=False, flag_value='mm', default='hs', show_default=True, help='Species.')
@click.option('-n', '--ntop', type=click.INT, default=10, show_default=True, help='Top terms to plot.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True))
def main(outdir, bname, condaenv, numthreads, species, ntop, infile):
	"""
Perform a pathway enrichment analysis using `pathfindR::run_pathfindR()`.

INFILE is a DESeq2 result file. See `scrnah5ad2degcontrastpairbydeseq2`.

\b
Example:
  f=$(parentsearch.sh -d test_data/deseq2pathfindr snRNA_AC_cluster2_HAC10VsHAC26.txt.gz)
  bname=$(basename "$f" .txt.gz)
  outdir=$(mktemp -d -u)
  deseq2pathfindr -d "$outdir" -b "$bname" -- "$f"

\b
See also:
  Upstream:
    scrnah5ad2degcontrastpairbydeseq2
  Depends:
    R/pathfindR

\b
Date: 2023/04/17
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f'{absdir}/R/{scriptname}.R'
	Path(outdir).mkdir(parents=True, exist_ok=True)
	os.chdir(outdir)

	pathfindf=f"{bname}_pathfind.txt.gz"
	degresult=pd.read_csv(infile, header=0, sep='\t')
	degresult[['symbol', 'log2FoldChange', 'padj']].to_csv(pathfindf, sep='\t', index=False)

	if species=='hs':
		geneset=['KEGG', 'GO-BP', 'GO-CC', 'GO-MF']
		pinpath='Biogrid'
		convert2alias='T'
	elif species=='mm':
		geneset=['mmu_KEGG']
		pinpath=['mmu_STRING']
		convert2alias='F'

	def cmd(*args):
		gset, =args
		exprs=[
			f"gset='{gset}'",
			f"bname='{bname}_{gset}'",
			f"pinpath='{pinpath}'",
			f"symbolalias={convert2alias}",
			f"ntop={ntop}",
			f"numthreads={numthreads}",
			f"infile='{pathfindf}'",
			]
		exeR.callback(exprs, script=script, condaenv=condaenv, verbose=True)

	with Pool(processes=numthreads) as p:
		p.starmap(cmd, zip(geneset))

	return 0

if __name__ == "__main__":
	sys.exit(main())
