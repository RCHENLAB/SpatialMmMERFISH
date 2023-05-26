#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import scanpy as sc
if f.endswith('.h5ad'):
	x=sc.read(f)
else:
	print("Error: please input the imputed .h5ad file.")
	import sys
	sys.exit(-1)

obs=x.obs
if len(obss)>0:
	obs=x.obs[obss]
print(obs)

import numpy as np
import pandas as pd
import anndata as ad
x=ad.AnnData(X=np.rint(x.X), obs=obs, var=pd.DataFrame(index=x.var.index))
sc.pp.filter_genes(x, min_counts=mincount)

sc.write(filename=f'{outdir}/{bname}.h5ad', adata=x)
x.obs['barcode']=x.obs.index
x.obs.to_csv(f'{outdir}/{bname}_obs.txt.gz', sep='\t', index=False)
x.var['symbol']=x.var.index
x.var.to_csv(f'{outdir}/{bname}_var.txt.gz', sep='\t', index=False)
