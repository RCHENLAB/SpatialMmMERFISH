#!/usr/bin/env python3
# vim: set noexpandtab tabstop=2 shiftwidth=2 softtabstop=-1 fileencoding=utf-8:

import pandas as pd
x=pd.read_csv(f, index_col=0)

if norm:
	x=x.div(x.sum(axis=1), axis=0)

from matplotlib import pyplot as plt
plt.figure(figsize=(width, height), dpi=500)
plt.hist(x.max(axis=1), bins=100)
plt.savefig(f'{outdir}/{bname}.pdf', bbox_inches='tight')
plt.close()
