#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

if __name__ == '__main__':
	import scanpy as sc
	import os
	if os.path.isdir(f):
		x=sc.read_10x_mtx(f)
	elif f.endswith('.h5'):
		x=sc.read_10x_h5(f)
	elif f.endswith('.h5ad'):
		x=sc.read(f)

	import pandas as pd
	feature=pd.read_csv(feature, sep='\t', header=None)

	feature=[g for g in feature[0] if g.lower() in x.var.index]
	symbol=[g.lower() for g in feature]
	x=x[:, symbol].copy()
	x.var.index=feature
	x.obs_names_make_unique()
	x.var_names_make_unique()

	sc.write(filename=f'{outdir}/{bname}.h5ad', adata=x)
	x.obs['barcode']=x.obs.index
	x.obs.to_csv(f'{outdir}/{bname}_obs.txt.gz', sep='\t', index=False)
	x.var['symbol']=x.var.index
	x.var.to_csv(f'{outdir}/{bname}_var.txt.gz', sep='\t', index=False)
