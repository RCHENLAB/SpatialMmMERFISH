import yaml
import subprocess
from pathlib import Path

# Parameters
infile, outdir, bname, nowtimestr=config['infile'].split(','), config['outdir'], config['bname'].split(','), config['nowtimestr']
config['infile'], config['bname']=infile, bname
files=dict(zip(bname, infile))

# for rules using a function
def get_file(wildcards):
	return files[wildcards.name]

if config['condaenv']=='None':
	config['condaenv']=None

# Default parameters
if 'deseq2enhancedvolcano' not in config:
	config['deseq2enhancedvolcano']={
		'height': 6,
		'width': 6,
		'qthr': 0.05,
		'log2fcthr': 1.0,
		'xlim': 10,
		'ylim': 15,
	}
deseq2enhancedvolcano_cmd=' '.join([
	f"deseq2enhancedvolcano.sh",
	f"-H {config['deseq2enhancedvolcano']['height']}",
	f"-W {config['deseq2enhancedvolcano']['width']}",
	f"-q {config['deseq2enhancedvolcano']['qthr']}",
	f"-l {config['deseq2enhancedvolcano']['log2fcthr']}",
	f"--xlim {config['deseq2enhancedvolcano']['xlim']}",
	f"--ylim {config['deseq2enhancedvolcano']['ylim']}",
	])



if 'rplothistogram' not in config:
	config['rplothistogram']={
		'height': 6,
		'width': 6,
		'breaks': 100,
		'xlab': 'P-value',
	}
rplothistogram_cmd=' '.join([
	f"rplothistogram.sh",
	f"-H {config['rplothistogram']['height']}",
	f"-W {config['rplothistogram']['width']}",
	f"-x {config['rplothistogram']['xlab']}",
	f"-b {config['rplothistogram']['breaks']}",
	])



if 'deseq2enrichgo' not in config:
	config['deseq2enrichgo']={
		'qthr': 0.05,
		'log2fcthr': 1.0,
		'species': 'hs',
	}
deseq2enrichgo_cmd=' '.join([
	f"deseq2enrichgo.sh",
	f"-q {config['deseq2enrichgo']['qthr']}",
	f"-l {config['deseq2enrichgo']['log2fcthr']}",
	f"-s {config['deseq2enrichgo']['species']}",
	])


if 'scrnaenrichgo2dotplot' not in config:
	config['scrnaenrichgo2dotplot']={
		'condaenv': config['condaenv'],
		'height': 5,
		'width': 6,
		'ntop': 10,
		'fontsize': 10,
	}
if config['condaenv']:
	config['scrnaenrichgo2dotplot']['condaenv']=config['condaenv']

scrnaenrichgo2dotplot_cmd=[
	f"scrnaenrichgo2dotplot",
	f"-H {config['scrnaenrichgo2dotplot']['height']}",
	f"-W {config['scrnaenrichgo2dotplot']['width']}",
	f"-n {config['scrnaenrichgo2dotplot']['ntop']}",
	f"-f {config['scrnaenrichgo2dotplot']['fontsize']}",
	]
if config['scrnaenrichgo2dotplot']['condaenv']:
	scrnaenrichgo2dotplot_cmd+=[f"-e {config['scrnaenrichgo2dotplot']['condaenv']}"]
scrnaenrichgo2dotplot_cmd=' '.join(scrnaenrichgo2dotplot_cmd)


if 'scrnaenrichgo2enrichplot' not in config:
	config['scrnaenrichgo2enrichplot']={
		'height': 6,
		'width': 10,
		'ntop': 30,
	}
scrnaenrichgo2enrichplot_cmd=' '.join([
	f"scrnaenrichgo2enrichplot.sh",
	f"-H {config['scrnaenrichgo2enrichplot']['height']}",
	f"-W {config['scrnaenrichgo2enrichplot']['width']}",
	f"-n {config['scrnaenrichgo2enrichplot']['ntop']}",
	])


if 'deseq2pathfindr' not in config:
	config['deseq2pathfindr']={
		'skip': True,
		'condaenv': config['condaenv'],
		'species': 'hs',
		'ntop': 10,
	}
if config['condaenv']:
	config['deseq2pathfindr']['condaenv']=config['condaenv']

deseq2pathfindr_cmd=[
	f"deseq2pathfindr",
	f"-s {config['deseq2pathfindr']['species']}",
	f"-n {config['deseq2pathfindr']['ntop']}",
	]
if config['deseq2pathfindr']['condaenv']:
	deseq2pathfindr_cmd+=[f"-e {config['deseq2pathfindr']['condaenv']}"]
deseq2pathfindr_cmd=' '.join(deseq2pathfindr_cmd)


# debug parameters
with open(f"config_{nowtimestr}.yaml", 'w') as f:
	yaml.dump(config, f, sort_keys=False)
