# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(rhdf5))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(DropletUtils))
suppressPackageStartupMessages(library(Matrix))

X=h5read(f, 'X')
print('==> X')
str(X)
## array to vector
obsindex=h5read(f, '/obs/_index', drop=T)
print('==> obs/_index')
str(obsindex)
varindex=h5read(f, '/var/_index', drop=T)
print('==> var/_index')
str(varindex)
colnames(X)=obsindex
rownames(X)=varindex
print('==> X')
str(X)

dds=DESeqDataSetFromMatrix(
	countData=X
	, colData=data.frame(obsindex=obsindex)
	, design=~1
	)
print('==> DESeqDataSetFromMatrix()')
str(dds)

cntraw=counts(dds, normalized=F)
print('==> cntraw')
str(cntraw)
write10xCounts(sprintf('%s/%s_raw.h5', outdir, bname), Matrix(cntraw, sparse=T))

pdf(sprintf('%s/%s_raw_hist_rowMeans.pdf', outdir, bname), width=width, height=height)
hist(rowMeans(cntraw), breaks=100, xlab='rowMeans of raw counts', main='')
dev.off()

pdf(sprintf('%s/%s_raw_hist_rowMeans_asinh.pdf', outdir, bname), width=width, height=height)
hist(asinh(rowMeans(cntraw)), breaks=100, xlab='asinh(rowMeans of raw counts)', main='')
dev.off()

# filtering genes by base mean without normalization
# may be useful to filter out extremely lowly expressed genes
if (rowmean>0) {
	dds=subset(dds, rowMeans(counts(dds, normalized=F))>=rowmean)
	dds=subset(dds, select=colSums(counts(dds, normalized=F))>0) ## excluding empty cells
}

tmp=computeSumFactors(dds, BPPARAM=MulticoreParam(numthreads))
sizeFactors(dds)=sizeFactors(tmp)
rm(tmp)
gc()
print('==> scran::computeSumFactors()')
str(dds)

pdf(sprintf('%s/%s_sizefactor_hist.pdf', outdir, bname), width=width, height=height)
hist(sizeFactors(dds), breaks=100)
dev.off()
write.txt(cbind(barcode=rownames(colData(dds)), colData(dds)), file=sprintf('%s/%s_colData.txt.gz', outdir, bname))


cntnorm=counts(dds, normalized=T)
print('==> cntnorm')
str(cntnorm)
write10xCounts(sprintf('%s/%s_norm.h5', outdir, bname), Matrix(cntnorm, sparse=T))

pdf(sprintf('%s/%s_norm_hist_rowMeans.pdf', outdir, bname), width=width, height=height)
hist(rowMeans(cntnorm), breaks=100, xlab='rowMeans of normalized counts', main='')
dev.off()

pdf(sprintf('%s/%s_norm_hist_rowMeans_asinh.pdf', outdir, bname), width=width, height=height)
hist(asinh(rowMeans(cntnorm)), breaks=100, xlab='asinh(rowMeans of normalized counts)', main='')
dev.off()
