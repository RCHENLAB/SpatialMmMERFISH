#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import pandas as pd
rnatype=pd.read_csv(rnatype, sep='\t', header=None)
merfishtype=pd.read_csv(merfishtype, sep='\t', header=None)
print(rnatype)
print(merfishtype)
rnatype=rnatype[0]
merfishtype=merfishtype[0]

import numpy as np
mask=np.tile(rnatype, (merfishtype.size, 1)).T==np.tile(merfishtype, (rnatype.size, 1))
print(mask)

tmp=np.sum(mask, axis=0)
if np.min(tmp)==0:
	print(f'Warning: empty MERFISH cells detected.')

import pandas as pd
x_freq=pd.DataFrame(np.asarray(np.unique(tmp, return_counts=True)).T, columns=['value', 'count'])
print(x_freq)
print(type(x_freq))
x_freq.to_csv(f'{outdir}/{bname}_freq.txt.gz', sep='\t', index=False, header=True)

from matplotlib import pyplot as plt
plt.figure(figsize=(5, 5), dpi=500)
plt.hist(tmp, bins=50)
plt.xlabel('Number of RNA cells per MERFISH cell')
plt.ylabel('Frequency')
plt.tight_layout()
plt.savefig(f'{outdir}/{bname}_hist.pdf', bbox_inches='tight')
plt.close()

import scanpy as sc
map=sc.read(map)
print(map.X)
map.X[~mask]=0
print(map.X)

if totalweight>0:
	map.X=totalweight*map.X/np.sum(map.X, axis=0)
	print(np.sum(map.X, axis=0))

def saveh5ad(adata, outdir, oname):
	if '_index' in adata.var_keys():
		adata.var.set_index('_index', inplace=True)
	if adata.raw is not None:
		adata.__dict__['_raw'].__dict__['_var']=adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
	sc.write(filename=f'{outdir}/{oname}.h5ad', adata=adata)
	adata.obs['barcode']=adata.obs.index
	adata.obs.to_csv(f'{outdir}/{oname}_obs.txt.gz', sep='\t', index=False)
	adata.var['symbol']=adata.var.index
	adata.var.to_csv(f'{outdir}/{oname}_var.txt.gz', sep='\t', index=False)

saveh5ad(map, outdir, f'{bname}_mapbymask')

rna=sc.read(rna)
print(rna)
sc.pp.filter_genes(rna, min_cells=1)
rna.X=rna.X.toarray()

X_space=map.X.T @ rna.X

print(X_space)
print(type(X_space))
print(X_space.shape)

import anndata as ad
x=ad.AnnData(X=X_space, obs=map.var, var=rna.var, uns=rna.uns)
sc.write(filename=f'{outdir}/{bname}.h5ad', adata=x)
x.obs['barcode']=x.obs.index
x.obs.to_csv(f'{outdir}/{bname}_obs.txt.gz', sep='\t', index=False)
x.var['symbol']=x.var.index
x.var.to_csv(f'{outdir}/{bname}_var.txt.gz', sep='\t', index=False)
