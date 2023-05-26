#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import pandas as pd
x=pd.read_csv(f, index_col=0)
xmax=x.max(axis=1)
x['celltype']=x.idxmax(axis=1)
x['barcode']=x.index
x['max']=xmax
print(x)
x[['barcode', 'max', 'celltype']].to_csv(f'{outdir}/{bname}.txt.gz', sep='\t', index=False)
