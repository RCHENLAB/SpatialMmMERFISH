#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import scanpy as sc
x=sc.read(f)
df=x.uns["train_genes_df"]

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns

fig, axs=plt.subplots(1, 4, figsize=(12, 3), sharey=True)
axs_f=axs.flatten()

axs_f[0].set_ylim([0.0, 1.0])
for i in range(1, len(axs_f)):
	axs_f[i].set_xlim([0.0, 1.0])
	axs_f[i].set_ylim([0.0, 1.0])
sns.histplot(data=df, y="train_score", bins=bins, ax=axs_f[0], color="coral")

axs_f[1].set_title("score vs sparsity (single cells)")
sns.scatterplot(
	data=df
	, y="train_score"
	, x="sparsity_sc"
	, ax=axs_f[1]
	, alpha=0.5
	, color="coral"
)

axs_f[2].set_title("score vs sparsity (spatial)")
sns.scatterplot(
	data=df
	, y="train_score"
	, x="sparsity_sp"
	, ax=axs_f[2]
	, alpha=0.5
	, color="coral"
)

axs_f[3].set_title("score vs sparsity (sp - sc)")
sns.scatterplot(
	data=df
	, y="train_score"
	, x="sparsity_diff"
	, ax=axs_f[3]
	, alpha=0.5
	, color="coral"
)

plt.tight_layout()
plt.savefig(f'{outdir}/{bname}.pdf')
plt.close()

df['symbol']=df.index
df.to_csv(f'{outdir}/{bname}_train_genes_df.txt.gz', sep='\t', index=False)

