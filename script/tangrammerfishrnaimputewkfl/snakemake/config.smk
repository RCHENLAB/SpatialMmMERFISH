import yaml
import json # for debug
import re

# Parameters
infile, outdir, bname, reference, nowtimestr=config['infile'].split(','), config['outdir'], config['bname'].split(','), config['reference'], config['nowtimestr']
config['infile'], config['bname']=infile, bname
name2file=dict(zip(bname, infile))

# for rules using a function
def get_file(wildcards):
	return name2file[wildcards.name]

# Default parameters
default_params={
	'tangramrnamerfish2integrate': {
		'condaenv': None,
		'rna': reference,
		'clusterlabel': 'majorclass',
		'seed': 12345,
		'gpu': [0, 1],
		'mode': 'cells',
		'cpu': False,
	},
	'tangramcellprojection2histogram': {
		'condaenv': None,
		'height': 5,
		'width': 5,
		'norm': False,
	},
	'tangrammap2plottrainingscores': {
		'condaenv': None,
		'bins': 50,
	},
	'tangrammerfishcelltypepred2label': {
		'condaenv': None,
	},
	'tangramintegrate2projectgenesfilterbycelltype': {
		'condaenv': None,
		'rnalabel': 'majorclass',
		'merfishlabel': 'bindSC_subtype',
		'totalweight': -1,
	},
	'tangramimpute2roundh5ad': {
		'condaenv': None,
		'mincount': 1,
	},
}

# Assign default params to config if not existing
for key, paramdict in default_params.items():
	if key not in config:
		config[key]=paramdict
	else:
		for subkey, value in paramdict:
			if subkey not in config[key]:
				config[key][subkey]=value
			else:
				print("Info: {key} and {subkey} are in config. So, skipping default assignment.")

# Overwrite condaenv
if config['condaenv']=='None':
	config['condaenv']=None
if config['condaenv']:
	for key, paramdict in config.items():
		if isinstance(paramdict, dict) and 'condaenv' in paramdict:
			paramdict['condaenv']=config['condaenv']

# Overwrite configurations by KEY=VALUE definition pairs
if config['define']:
	for kv in config['define'].split(','):
		m=re.match(r"^(\w+)\.?(\w+)?:([=sifbSIFB]):(\S+)$", kv)
		if m:
			key, subkey, op, value=m.groups()
			if subkey:
				config[key][subkey]=opvalue(op, value)
			else:
				config[key]=opvalue(op, value)
		else:
			print(f"Warning: {kv} not recognized.")
print(json.dumps(config, indent=4))

# save parameters
with open(f"config_{nowtimestr}.yaml", 'w') as f:
	yaml.dump(config, f, sort_keys=False)

# Do the work
# Build command string for each step
# +
tangramrnamerfish2integrate_cmd=[
	f"tangramrnamerfish2integrate.sh",
	f"-r {config['tangramrnamerfish2integrate']['rna']}",
	f"-l {config['tangramrnamerfish2integrate']['clusterlabel']}",
	f"-s {config['tangramrnamerfish2integrate']['seed']}",
	f"-m {config['tangramrnamerfish2integrate']['mode']}",
	]
if config['tangramrnamerfish2integrate']['condaenv']:
	tangramrnamerfish2integrate_cmd+=[f"-e {config['tangramrnamerfish2integrate']['condaenv']}"]
tangramrnamerfish2integrate_cmd=' '.join(tangramrnamerfish2integrate_cmd)
# -

# +
tangramcellprojection2histogram_cmd=[
	f"tangramcellprojection2histogram.sh",
	f"-H {config['tangramcellprojection2histogram']['height']}",
	f"-W {config['tangramcellprojection2histogram']['width']}",
	]
if config['tangramcellprojection2histogram']['condaenv']:
	tangramcellprojection2histogram_cmd+=[f"-e {config['tangramcellprojection2histogram']['condaenv']}"]
if config['tangramcellprojection2histogram']['norm']:
	tangramcellprojection2histogram_cmd+=[f"-n"]
tangramcellprojection2histogram_cmd=' '.join(tangramcellprojection2histogram_cmd)
# -

# +
tangrammap2plottrainingscores_cmd=[
	f"tangrammap2plottrainingscores.sh",
	f"-k {config['tangrammap2plottrainingscores']['bins']}",
	]
if config['tangrammap2plottrainingscores']['condaenv']:
	tangrammap2plottrainingscores_cmd+=[f"-e {config['tangrammap2plottrainingscores']['condaenv']}"]
tangrammap2plottrainingscores_cmd=' '.join(tangrammap2plottrainingscores_cmd)
# -

# +
tangrammerfishcelltypepred2label_cmd=[
	f"tangrammerfishcelltypepred2label.sh",
	]
if config['tangrammerfishcelltypepred2label']['condaenv']:
	tangrammerfishcelltypepred2label_cmd+=[f"-e {config['tangrammerfishcelltypepred2label']['condaenv']}"]
tangrammerfishcelltypepred2label_cmd=' '.join(tangrammerfishcelltypepred2label_cmd)
# -

# +
tangramintegrate2projectgenesfilterbycelltype_cmd=[
	f"tangramintegrate2projectgenesfilterbycelltype.sh",
	f"-w {config['tangramintegrate2projectgenesfilterbycelltype']['totalweight']}",
	]
if config['tangramintegrate2projectgenesfilterbycelltype']['condaenv']:
	tangramintegrate2projectgenesfilterbycelltype_cmd+=[f"-e {config['tangramintegrate2projectgenesfilterbycelltype']['condaenv']}"]
tangramintegrate2projectgenesfilterbycelltype_cmd=' '.join(tangramintegrate2projectgenesfilterbycelltype_cmd)
# -

# +
tangramimpute2roundh5ad_cmd=[
	f"tangramimpute2roundh5ad.sh",
	f"-m {config['tangramimpute2roundh5ad']['mincount']}",
	]
if config['tangramimpute2roundh5ad']['condaenv']:
	tangramimpute2roundh5ad_cmd+=[f"-e {config['tangramimpute2roundh5ad']['condaenv']}"]
tangramimpute2roundh5ad_cmd=' '.join(tangramimpute2roundh5ad_cmd)
# -
