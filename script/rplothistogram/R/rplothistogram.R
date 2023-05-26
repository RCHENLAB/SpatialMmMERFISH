# vim: set noexpandtab tabstop=2:

f=read.table('stdin', header=F, sep='\t', quote="", check.names=F, comment.char='#', stringsAsFactors=F)[,1]

if (transform) {
	f=asinh(f)
}

if (fixedwidth) {
	breaks=seq(lower, higher, by=(higher-lower)/breaks)
}

pdf(outfile, width=width, height=height)
if (combine) {
	x=hist(f, breaks=breaks, xlab=xlab, main=main, freq=!density, right=right)
	x=hist(f, breaks=breaks, xlab=xlab, main=main, freq=density, right=right)
} else {
	x=hist(f, breaks=breaks, xlab=xlab, main=main, freq=!density, right=right)
}
dev.off()
write.table(data.frame(mids=x$mids, density=x$density, counts=x$counts), file=stdout(), sep='\t', quote=F, row.names=F, col.names=T)
