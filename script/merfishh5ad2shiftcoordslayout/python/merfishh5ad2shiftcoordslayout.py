#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import scanpy as sc
import pandas as pd
import numpy as np

x=sc.read(f)
print(x)

keys=x.obs[shiftby].unique()
if len(keys)<2:
	print(keys)
	print(f'{bname}: skipping shift due to too few groups.')
	import sys
	sys.exit(-1)

# centerize
def centerize(obs):
	tmp=obs[[coordx, coordy]]
	tmp-=np.amin(tmp, axis=0)
	return tmp
x.obs[['centerize_x', 'centerize_y']]=x.obs.groupby(shiftby).apply(centerize)
print(x.obs.describe())

codes, _=pd.factorize(x.obs[shiftby])
shiftmatrix=np.multiply(np.column_stack((codes//ncol, codes%ncol)), [stepx, stepy])
print(shiftmatrix)

x.obs[[finalx, finaly]]=np.add(x.obs[['centerize_x', 'centerize_y']], shiftmatrix)
print(x)

sc.write(filename=f'{outdir}/{bname}.h5ad', adata=x)
x.obs['barcode']=x.obs.index
x.obs.to_csv(f'{outdir}/{bname}_obs.txt.gz', sep='\t', index=False)
x.var['symbol']=x.var.index
x.var.to_csv(f'{outdir}/{bname}_var.txt.gz', sep='\t', index=False)
