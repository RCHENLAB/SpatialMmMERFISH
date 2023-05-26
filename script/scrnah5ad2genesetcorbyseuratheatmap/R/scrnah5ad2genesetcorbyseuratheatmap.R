# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(pheatmap))

x=h5ad2seurat(infile)

if (normonly) {
	x=x %>% NormalizeData()
	X=as.matrix(x[['RNA']]@data)
	print('==> X, normonly')
	str(X)
} else {
	x=x %>% NormalizeData() %>% ScaleData(features=rownames(x))
	X=x[['RNA']]@scale.data
	print('==> X, scaled')
	str(X)
}

genes=read.txt(genefile, header=F)[, 1]
print('==> genes')
str(genes)

genes=genes[genes%in%rownames(X)]
print('==> genes, in features')
str(genes)

if (is.null(vfilter)) {
	dm=X[genes, , drop=F]
	PCC=cor(t(dm), method=method)
} else {
	PCC=matrix(NA, length(genes), length(genes))
	rownames(PCC)=genes
	colnames(PCC)=genes
	lapply(
		combn(genes, 2, simplify=F)
		, function(twogenes) {
			g1=twogenes[1]
			g2=twogenes[2]
			dm=t(X[twogenes, ])
			dm=subset(dm, (dm[, g1]>=vfilter & dm[, g2]>=vfilter))
			pcc=cor(dm[, g1], dm[, g2], method=method)
			PCC[g1, g2]<<-pcc
			PCC[g2, g1]<<-pcc
		})
	diag(PCC)=1
}
print('==> PCC')
print(PCC)
str(PCC)
write.txt(PCC, file=sprintf('%s/%s.txt.gz', outdir, bname), row.names=T, col.names=T)

if (scatter) {
	apply(
		expand.grid(rownames(PCC), colnames(PCC))
		, 1
		, function(twogenes) {
			g1=twogenes[1]
			g2=twogenes[2]
			dm=t(X[twogenes, ])

			if (!is.null(vfilter)) {
				dm=subset(dm, (dm[, g1]>=vfilter & dm[, g2]>=vfilter))
			}

			pdf(sprintf('%s/%s_%s_%s.pdf', outdir, bname, g1, g2), width=5, height=5)
			smoothScatter(
				dm[, g1]
				, dm[, g2]
				, colramp = colorRampPalette(c('blue', 'yellow', 'red'))
				, xlab=g1
				, ylab=g2
				, main=sprintf('%s Vs %s: PCC=%.2f', g1, g2, PCC[g1, g2])
				, asp=NULL
				)
			dev.off()
		})
}

pdf(sprintf('%s/%s.pdf', outdir, bname), width=width, height=height)
x=pheatmap(
	PCC
	, color=colorRampPalette(c('blue', 'white', 'red'))(1024)
	, breaks=seq(-colorkey, colorkey, length.out=1024)
	, legend_breaks=c(-colorkey, 0, colorkey)
	, legend_labels=c(-colorkey, 0, colorkey)
	, show_rownames=T
	, show_colnames=T
	, cluster_rows=T
	, cluster_cols=T
	, display_numbers=shownumber
	)
dev.off()
