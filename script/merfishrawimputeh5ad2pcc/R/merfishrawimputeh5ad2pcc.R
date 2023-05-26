# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(parallel))

merfishraw=h5ad2seurat(rawfile)
merfishraw=merfishraw[['RNA']]@counts
print('==> merfishraw')
str(merfishraw)

merfishimpute=h5ad2seurat(imputefile)
merfishimpute=merfishimpute[['RNA']]@counts
print('==> merfishimpute')
str(merfishimpute)

commonsymbol=intersect(rownames(merfishraw), rownames(merfishimpute))
commoncell=intersect(colnames(merfishraw), colnames(merfishimpute))

if (feature!='') {
	feature=read.txt(feature, header=F)[, 1]
	commonsymbol=intersect(commonsymbol, feature)
}

if (barcode!='') {
	barcode=read.txt(barcode, header=F)[, 1]
	commoncell=intersect(commoncell, barcode)
}

print('========> commonsymbol')
str(commonsymbol)
print('========> commoncell')
str(commoncell)

merfishraw=merfishraw[commonsymbol, commoncell]
print('========> merfishraw, common')
str(merfishraw)

merfishimpute=merfishimpute[commonsymbol, commoncell]
print('========> merfishimpute, common')
str(merfishimpute)

result=do.call(rbind
	, parallel::mclapply(
		commonsymbol
		, function (symbol) {
			x=merfishimpute[symbol, ]
			y=merfishraw[symbol, ]
			pcc=cor(x, y, method='pearson')
			data.frame(symbol=symbol, pcc=pcc)
		}
		, mc.cores=numthreads
		)
	)
str(result)

write.txt(result, file=sprintf('%s/%s.txt.gz', outdir, bname))
