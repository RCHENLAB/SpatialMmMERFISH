# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(pheatmap))
x=read.txt(f, header=T)
print('==> x')
str(x)

if (label!='') {
	label=read.txt(label, header=F)[, 1]
	print('==> label')
	str(label)

	label=label[label %in% rownames(x)]
	print('==> label, common')
	str(label)

	x=x[label, label]
	print('==> x, reorder')
	str(x)
}

pdf(sprintf('%s/%s.pdf', outdir, bname), width=width, height=height)
x=pheatmap(
	x
	, color=colorRampPalette(c('blue', 'white', 'red'))(1024)
	, breaks=seq(-colorkey, colorkey, length.out=1024)
	, legend_breaks=c(-colorkey, 0, colorkey)
	, legend_labels=c(-colorkey, 0, colorkey)
	, show_rownames=T
	, show_colnames=T
	, cluster_rows=cluster
	, cluster_cols=cluster
	, display_numbers=shownumber
	)
dev.off()
