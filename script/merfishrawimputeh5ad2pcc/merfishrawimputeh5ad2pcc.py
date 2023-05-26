#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import sys
from exeR.exeR import exeR
from pathlib import Path
import click
CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(), default='.', show_default=True, help='Outdir.')
@click.option('-b', '--bname', type=click.STRING, required=True, help='Bname.')
@click.option('-r', '--rawfile', type=click.Path(exists=True), required=True, help='MERFISH raw count file.')
@click.option('-f', '--feature', type=click.Path(exists=True), help='Subset of feature in a file.')
@click.option('-c', '--barcode', type=click.Path(exists=True), help='Subset of barcode in a file.')
@click.option('-t', '--numthreads', type=click.INT, default=8, show_default=True, help='Number of threads.')
@click.argument('infile', type=click.Path(exists=True))
def main(outdir, bname, rawfile, feature, barcode, numthreads, infile):
	"""
To calculate Pearson Correlation Coefficient between imputation and MERFISH raw counts.

INFILE is a .h5ad file for imputation. See tangrammerfishrnaimputewkfl.

\b
Example:
  raw=$(parentsearch.sh -d test_data/merfishrawimputeh5ad2pcc merfish_raw.h5ad)
  impute=$(parentsearch.sh -d test_data/merfishrawimputeh5ad2pcc merfish_impute.h5ad)
  feature=$(parentsearch.sh -d test_data/merfishrawimputeh5ad2pcc Rod.txt.gz)
  barcode=$(parentsearch.sh -d test_data/merfishrawimputeh5ad2pcc VZG105a_WT2_Rod.txt.gz)
  bname=merfish_pcc
  outdir=$(mktemp -d -u)
  merfishrawimputeh5ad2pcc -d "$outdir" -b "$bname" -r "$raw" -f "$feature" -c "$barcode" -- "$impute"

\b
See also:
  Upstream:
    tangrammerfishrnaimputewkfl
  Depends:
    R/Seurat

\b
Date: 2023/05/23
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	Path(outdir).mkdir(parents=True, exist_ok=True)
	absdir=Path(__file__).parent
	scriptname=Path(__file__).stem
	script=f'{absdir}/R/{scriptname}.R'
	exprs=[
		f"outdir='{outdir}'",
		f"bname='{bname}'",
		f"rawfile='{rawfile}'",
		f"imputefile='{infile}'",
		f"feature='{feature if feature is not None else ''}'",
		f"barcode='{barcode if barcode is not None else ''}'",
		f"numthreads={numthreads}",
		]
	return exeR.callback(exprs, script=script, condaenv=None, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
