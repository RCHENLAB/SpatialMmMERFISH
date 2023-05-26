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
@click.option('-D', '--define', type=click.STRING, multiple=True, help='Define KEY=VALUE pairs.')
@click.option('-r', '--reference', type=click.Path(exists=False, resolve_path=True), help='Reference single cell RNA .h5ad file.')
@click.option('-t', '--numthreads', type=click.INT, default=4, show_default=True, help='Number of threads.')
@click.option('-b', '--bname', multiple=True, type=click.STRING, help='Bnames for infile. Default: basenames of infiles.')
@click.option('-n', '--dryrun', is_flag=True, help='Dry-run.')
@click.argument('infile', type=click.Path(exists=True, resolve_path=True), required=True, nargs=-1)
def main(outdir, condaenv, configfile, define, reference, numthreads, bname, dryrun, infile):
	"""
To perform imputation for MERFISH from RNA by Tangram.

INFILE are MERFISH .h5ad files.

\b
Example:
  files=(
  /storage/singlecell/jinli/wkfl/spatialprj/VA45/imputation/preproc/scrnah5adsplitby/VZG105a_VA45_merged_HC_VZG105a_WT1.h5ad
  /storage/singlecell/jinli/wkfl/spatialprj/VA45/imputation/preproc/scrnah5adsplitby/VZG105a_VA45_merged_HC_VZG105a_WT2.h5ad
  )
  ref=/storage/singlecell/jinli/wkfl/spatialprj/VA45/imputation/ref/download/HC.h5ad
  outdir=$(mrrdir.sh)
  # slurmtaco.sh -n g01 -t 2 -m 60G -- tangrammerfishrnaimputewkfl -d "$outdir" -e u_tangram -t 4 -r "$ref" -n -- "${files[@]}"
  slurmtaco.sh -n g01 -t 2 -m 60G -- tangrammerfishrnaimputewkfl -d "$outdir" -e u_tangram -t 4 -r "$ref" -c config.yaml -- "${files[@]}"

\b
  # A running example using `-D|--define`
  merfishdir=/storage/singlecell/jinli/wkfl/spatialprj/VA45/imputation/displacedAC/merfish/scrnah5adfiltergenescellsbycounts
  rnadir=/storage/singlecell/jinli/wkfl/spatialprj/VA45/imputation/displacedAC/rna
  outdir=$(mrrdir.sh)
  samples=(DV1 DV2 DV3 TN1 TN2 TN3 VZG105a_WT1 VZG105a_WT2 VZG105a_WT3 VZG105a_WT4)
  celltypes=(AC_17 AC_2 AC_7 AC_21 AC_44 AC_48 AC_32 AC_39 AC_36 AC_37 AC_54 AC_46)
  function cmd {
  local sampleid=$1
  local celltype=$2
  local merfish=$merfishdir/${sampleid}_${celltype}.h5ad
  local rna=$rnadir/${celltype}.h5ad
  local bname=$(basename "$merfish" .h5ad)
  local gpu=$((RANDOM%4))
  local define=(
  tangramrnamerfish2integrate.clusterlabel:=:subclass_label
  tangramrnamerfish2integrate.gpu:I:"$gpu"
  tangramintegrate2projectgenesfilterbycelltype.rnalabel:=:subclass_label
  tangramintegrate2projectgenesfilterbycelltype.merfishlabel:=:bindSC_subtype
  )
  if fileexists.sh "$merfish" "$rna"
  then
  	slurmtaco.sh -j "$bname" -b -- tangrammerfishrnaimputewkfl -d "$outdir/$bname" -e u_tangram -t 2 -r "$rna" $(basharr2cmdopts.sh -o -D -- "${define[@]}") -- "$merfish"
  fi
  }
  source env_parallel.bash
  env_parallel -j 8 cmd ::: "${samples[@]}" ::: "${celltypes[@]}"

\b
Note:
  1. `-D|--define` for KEY=VALUE pairs.

\b
  ```
  m=re.match(r"^(\w+)\.?(\w+)?:([=sifbSIFB]):(\S+)$", kv)
  ```
  Format: key[.subkey]:operation:value

\b
  1.1 Operations
  - `=` or `s`: a direct string assignment
  - `S`: A string list assignment. The elements are separated by `::`.
  - `i`, `f`, `b`: type convertion using int(), float(), and bool(), respectively
  - `I`, `F`, `B`: List separated by `::`, and each element is converted by int(), float(), and bool(), respectively

\b
See also:
  Steps:
    tangramrnamerfish2integrate.sh
    tangramcellprojection2histogram.sh
    tangrammap2plottrainingscores.sh
    tangrammerfishcelltypepred2label.sh
    tangramintegrate2projectgenesfilterbycelltype.sh
    tangramimpute2roundh5ad.sh
  Related:
    tangramtmpl.sh

\b
Date: 2023/05/22
Authors: Jin Li <lijin.abc@gmail.com>
	"""
	if len(bname)!=len(infile):
		bname=[Path(file).stem for file in infile]
	Path(outdir).mkdir(parents=True, exist_ok=True)
	absdir=Path(__file__).parent
	nowtimestr=datetime.datetime.now().strftime('%y%m%d_%H%M%S')

	for file in [
			f"{absdir}/snakemake/Snakefile",
			f"{absdir}/snakemake/def.smk",
			f"{absdir}/snakemake/config.smk",
			]:
		shutil.copy2(file, outdir) # copy Snakfile and config.smk

	# common command string
	cmdstr=[
		f"snakemake",
		f"-C outdir='{outdir}' reference='{reference}' infile='{','.join(infile)}' bname='{','.join(bname)}' condaenv='{condaenv}' nowtimestr='{nowtimestr}' define='{','.join(define)}'",
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
