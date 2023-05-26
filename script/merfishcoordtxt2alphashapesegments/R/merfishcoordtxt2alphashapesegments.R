# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(alphahull))
x=read.txt(f)
print('==> x')
str(x)

xcoord=x[, c(coordx, coordy)]
rownames(xcoord)=x[, 'barcode']
print('==> xcoord')
str(xcoord)
print(head(xcoord))

alphashape=ashape(as.matrix(xcoord), alpha=alpha)
print('==> alpha shape')
str(alphashape)
saveRDS(alphashape, file=sprintf('%s/%s.rds', outdir, bname))

extpts=alphashape$alpha.extremes
print('==> extreme points')
str(extpts)
print(head(extpts))
write.txt(
	cbind(
		index=extpts
		, barcode=rownames(xcoord[extpts, ])
		, xcoord[extpts, ]
		)
	, file=sprintf('%s/%s_extremepts.txt.gz', outdir, bname)
	)

edges=alphashape$edges
print('==> edges')
str(edges)
print(head(edges))
write.txt(edges, file=sprintf('%s/%s_segments.txt.gz', outdir, bname))

pdf(sprintf('%s/%s.pdf', outdir, bname), width=width, height=height)
# plot(alphashape, col=c(4, 1), xlab='x-coordinate', ylab='y-coordinate', pch='.', main=sprintf('alpha: %s', alpha), wpoints=T)
plot(xcoord, pch='.', main=sprintf('alpha: %s', alpha))
for (i in 1:nrow(edges)) {
	# use index for validation
	# segments(edges[i, 'x1'], edges[i, 'y1'], edges[i, 'x2'], edges[i, 'y2'], col=6)
	p1=xcoord[edges[i, 'ind1'], ]
	p2=xcoord[edges[i, 'ind2'], ]
	lines(rbind(p1, p2), col=6)
}

points(xcoord[extpts,], col=5, pch='.')

xoff=15
yoff=0
text(xcoord[extpts, 1]+xoff, xcoord[extpts, 2]+yoff, extpts, col=4, cex=0.1)
dev.off()
