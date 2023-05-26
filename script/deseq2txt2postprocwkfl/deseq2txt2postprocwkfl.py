#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import os
import sys
import shutil
from exeBash.exeBash import exeBash
from pathlib import Path
import datetime
import click

CONTEXT_SETTINGS=dict(help_option_names=['-h', '--help'])
@click.command(context_settings=CONTEXT_SETTINGS)
@click.option('-d', '--outdir', type=click.Path(resolve_path=True), default='.', show_default=True, help='Outdir.')
@click.option('-e', '--condaenv', type=click.STRING, help='Conda environment.')
@click.option('-c', '--configfile', type=click.Path(exists=False, resolve_path=True), help='Configuration file in the YAML format.')
@click.option('-t', '--numthreads', type=click.INT, default=4, show_default=True, help='Number of threads.')
@click.option('-b', '--bname', multiple=True, type=click.STRING, help='Bnames for infile. Default: basenames of infiles.')
@click.option('-n', '--dryrun', is_flag=True, help='Dry-run.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True), required=True, nargs=-1)
def main(outdir, condaenv, configfile, numthreads, bname, dryrun, infile):
	"""
To run postprocessing for DESeq2 TXT result files.

INFILE are .txt.gz files.

\b
Example:
  indir=$(mrrdir.sh ..)
  outdir=$(mrrdir.sh)
  shopt -s extglob
  slurmtaco.sh -t 2 -m 20G -n d00 -- deseq2txt2postprocwkfl -d "$outdir" -t 20 -n -- "$indir"/!(*colData).txt.gz

\b
See also:
  Upstream:
    scrnah5ad2bulkdegcondtypebydeseq2
    scrnah5ad2bulkdegtwofactorlrtbydeseq2
    scrnah5ad2degcontrastpairbydeseq2
  Steps:
    deseq2enhancedvolcano.sh
    rplothistogram.sh
    deseq2enrichgo.sh
    scrnaenrichgo2dotplot
    scrnaenrichgo2enrichplot.sh
    deseq2pathfindr
\b
Date: 2023/05/09
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	if len(bname)!=len(infile):
		# bname=[Path(file).stem for file in infile]
		bname=[Path(file).name.removesuffix(''.join(Path(file).suffixes)) for file in infile]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	absdir=Path(__file__).parent
	nowtimestr=datetime.datetime.now().strftime('%y%m%d_%H%M%S')

	for file in [
			f"{absdir}/snakemake/Snakefile",
			f"{absdir}/snakemake/config.smk",
			]:
		shutil.copy2(file, outdir) # copy Snakfile and config.smk

	# common command string
	cmdstr=[
		f"snakemake",
		f"-C outdir='{outdir}' infile='{','.join(infile)}' bname='{','.join(bname)}' condaenv='{condaenv}' nowtimestr='{nowtimestr}'",
		f"-d {outdir}",
		f"-j {numthreads}",
		]
	if configfile:
		cmdstr+=[f"--configfile {configfile}"]

	if dryrun: # dry-run
		cmdstr+=[f"-n -p"]

	else: # running
		cmdstr+=[
			f"-r -p --debug-dag",
			f"--stats Snakefile_{nowtimestr}.stats",
			]

	os.chdir(outdir)
	return exeBash.callback(cmdstr=cmdstr, condaenv=None, verbose=True)

if __name__ == "__main__":
	sys.exit(main())
