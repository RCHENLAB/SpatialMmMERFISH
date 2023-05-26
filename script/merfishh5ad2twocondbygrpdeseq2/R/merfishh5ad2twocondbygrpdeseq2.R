# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(BiocParallel))

dds=h5ad2deseq2(infile, variable=c(group, cond), design=as.formula(sprintf('~%s', cond)))
print('==> h5ad2deseq2')
str(dds)

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

metadata=as.data.frame(colData(dds))
parallel::mclapply(
	unique(metadata[, group])
	, function(g) {
		gdds=subset(dds, select=((metadata[, group]==g)&(metadata[, cond]%in%c(condA, condB))))
		if (ncol(gdds)>0) {
			colData(gdds)[, cond]=relevel(colData(gdds)[, cond], ref=condB)
			colData(gdds)[, cond]=droplevels(colData(gdds)[, cond])
			gdds=DESeq(
				gdds
				, fitType='glmGamPoi'
				, test='LRT'
				, useT=T
				, minmu=1e-6
				, minReplicatesForReplace=Inf
				, reduced=~1
				, parallel=F
				)
			cat('==> DESeq()', g, '\n')
			str(gdds)
			saveRDS(gdds, file=sprintf('%s/%s_%s_DESeq.rds', outdir, bname, g))
			res=data.frame(results(gdds))
			res=res[order(res$padj), ]
			write.txt(cbind(symbol=rownames(res), res), file=sprintf('%s/%s_%s.txt.gz', outdir, bname, g))
		} else {
			write(sprintf('Error: empty cells in the group %s.\n', g), stderr())
		}
	}
	, mc.cores=numthreads
	)
