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

extpts=chull(xcoord)
print('==> chull')
str(extpts)
print(head(extpts))
write.txt(
	cbind(
		index=extpts
		, barcode=rownames(xcoord[extpts, ])
		, xcoord[extpts, ])
	, file=sprintf('%s/%s_chull.txt.gz', outdir, bname)
	)

extpts=c(extpts, extpts[1])
curvesegs=do.call(
	rbind
	,	apply(
		cbind(extpts[-length(extpts)], extpts[-1])
		, 1
		, function(tmp) {
			p1=xcoord[tmp[1], ]
			p2=xcoord[tmp[2], ]
			distance=unname(sqrt((p1[1]-p2[1])^2 + (p1[2]-p2[2])^2))
			data.frame(idx1=tmp[1], idx2=tmp[2], distance=distance)
		})
	)
# remove the max distance, the longest segment is in the bottom of the Arc
curvesegs=curvesegs[-which.max(curvesegs[, 'distance']), , drop=F]
# save the segments
write.txt(curvesegs, file=sprintf('%s/%s_segments.txt.gz', outdir, bname))

curvesegs=apply(
	subset(curvesegs, select=-distance)
	, 1
	, function(tmp) {
		xcoord[tmp, ]
	})

pdf(sprintf('%s/%s.pdf', outdir, bname), width=width, height=height)
plot(xcoord, pch='.')
for (seg in curvesegs) {
	lines(seg, col=6)
}
points(xcoord[extpts,], col=5)

xoff=50
yoff=0
text(xcoord[extpts, 1]+xoff, xcoord[extpts, 2]+yoff, extpts, col=4, cex=0.3)
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
