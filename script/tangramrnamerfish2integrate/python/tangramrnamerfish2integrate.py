#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import scanpy as sc
if merfish.endswith('.h5'):
	merfish=sc.read_10x_h5(merfish)
elif merfish.endswith('.h5ad'):
	merfish=sc.read(merfish)
print(merfish)

rna=sc.read(rna)
rnaraw=rna
sc.pp.normalize_total(rna)
print(rna)
print(rna.obs[clusterlabel].value_counts())

import tangram as tg
tg.pp_adatas(rna, merfish)

cluster_label=None
if mode=='clusters':
	cluster_label=clusterlabel

target_count=None
if mode=='constrained':
	target_count=merfish.obs.shape[0]

x=tg.map_cells_to_space(
	rna
	, merfish
	, mode=mode
	, cluster_label=cluster_label
	, target_count=target_count
	, device=device
	, random_state=seed
	, density_prior='uniform'
	)
print(x)

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

saveh5ad(x, outdir, f'{bname}_map')
saveh5ad(rna, outdir, f'{bname}_rna')
saveh5ad(merfish, outdir, f'{bname}_merfish')

# merfish.obsm['tangram_ct_pred'] may not be saved correctly, so save it only
tg.project_cell_annotations(x, merfish, annotation=clusterlabel)
merfish.obsm['tangram_ct_pred'].to_csv(f'{outdir}/{bname}_merfish_tangram_ct_pred.csv', index=True)

xx=tg.project_genes(adata_map=x, adata_sc=rna, cluster_label=cluster_label)
saveh5ad(xx, outdir, f'{bname}_impute_norm')

xx=tg.project_genes(adata_map=x, adata_sc=rnaraw, cluster_label=cluster_label)
saveh5ad(xx, outdir, f'{bname}_impute_raw')
