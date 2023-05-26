# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(maptools))
x=read.txt(f)
print('==> x')
str(x)

xcoord=x[, c(coordx, coordy)]
rownames(xcoord)=x[, 'barcode']
print('==> xcoord')
str(xcoord)
print(head(xcoord))

curvesegs=read.txt(segment)
print('==> curvesegs')
str(curvesegs)

if (length(label)==4) { # column labels, (x1,y1,x2,y2)
	curvesegs=apply(
		subset(curvesegs, select=label)
		, 1
		, function(tmp) {
			matrix(tmp, 2, 2, byrow=T)
		}
		, simplify=F
		)
} else { # point indices
	curvesegs=apply(
		subset(curvesegs, select=1:2)
		, 1
		, function(tmp) {
			xcoord[tmp, ]
		})
}

print('==> curvesegs, coord')
str(curvesegs)

xlim=c(
	min(xcoord[, 1], sapply(curvesegs, function(seg) { min(seg[, 1]) }))
	, max(xcoord[, 1], sapply(curvesegs, function(seg) { max(seg[, 1]) }))
	)
ylim=c(
	min(xcoord[, 2], sapply(curvesegs, function(seg) { min(seg[, 2]) }))
	, max(xcoord[, 2], sapply(curvesegs, function(seg) { max(seg[, 2]) }))
	)

pdf(sprintf('%s/%s.pdf', outdir, bname), width=width, height=height)
plot(xcoord, pch='.', xlim=xlim, ylim=ylim)
for (seg in curvesegs) {
	lines(seg, col=6)
}
dev.off()

mindist2segments=function(p, segs) {
	tmp=do.call(
		rbind
		, lapply(
			segs
			, function(seg) {
				maptools::nearestPointOnSegment(seg, p)
			})
		)
	min(tmp[, 'distance'])
}

res=do.call(
	c
	, parallel::mclapply(
		seq_len(nrow(xcoord))
		, function(cellindex) {
			if (cellindex%%1000==1) {
				print(cellindex)
			}
			mindist2segments(unlist(xcoord[cellindex, ]), curvesegs)
		}
		, mc.cores=numthreads
		)
	)
print('==> res')
str(res)
write.txt(cbind(x, distance=res), file=sprintf('%s/%s.txt.gz', outdir, bname))
