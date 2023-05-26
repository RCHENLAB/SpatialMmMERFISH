#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import sys
import scanpy as sc
import os
import anndata as ad

def readfile(args):
	id=args[0]
	f=args[-1]
	print(f'start: {id}')
	if os.path.isdir(f):
		x=sc.read_10x_mtx(f)
	elif f.endswith('.h5'):
		x=sc.read_10x_h5(f)
	elif f.endswith('.h5ad'):
		x=sc.read(f)
	x.var_names_make_unique()
	if len(args)>2:
		for h, m in zip(header[1:-1], args[1:-1]):
			x.obs[h]=m
	print(f'end: {id}')
	return (id, x)

def readrawadata(args):
	id=args[0]
	f=args[-1]
	print(f'start: {id}')
	if f.endswith('.h5ad'):
		x=sc.read(f)
	else:
		print(f'Error: raw feature is not supported for {f}')
		sys.exit(-1)
	x.var_names_make_unique()

	# bug fix gene symbol using x.var
	x=ad.AnnData(X=x.raw.X, obs=x.obs, var=x.var)

	if len(args)>2:
		for h, m in zip(header[1:-1], args[1:-1]):
			x.obs[h]=m
	print(f'end: {id}')
	return (id, x)

if __name__ == '__main__':
	fnameinfo, outdir, bname, join, numthreads, raw, zerofill=sys.argv[1:]
	numthreads=int(numthreads)
	raw=raw=="True"
	fillvalue=0 if zerofill=="True" else None

	import pandas as pd
	from multiprocessing import Pool

	f=pd.read_table(fnameinfo, header=0)
	header=f.columns.tolist()
	print(f)
	if f.shape[1]>=2:
		f=list(f.itertuples(index=False, name=None))
	else:
		print("Error: please input at least two columns 'id[<TAB>metadata]<TAB>infile'.")
		sys.exit(-1)

	with Pool(processes=numthreads) as p:
		if raw:
			x=p.map(readrawadata, f)
		else:
			x=p.map(readfile, f)
	x=ad.concat(adatas=dict(x), axis=0, join=join, label=header[0], merge=None, fill_value=fillvalue)
	x.obs_names_make_unique()
	sc.write(filename=f'{outdir}/{bname}.h5ad', adata=x)
	x.obs['barcode']=x.obs.index
	x.obs.to_csv(f'{outdir}/{bname}_obs.txt.gz', sep='\t', index=False)
	x.var['symbol']=x.var.index
	x.var.to_csv(f'{outdir}/{bname}_var.txt.gz', sep='\t', index=False)
