# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(BiocParallel))

dds=h5ad2deseq2(infile, variable=c(group, cond), design=as.formula(sprintf('~%s+%s', group, cond)))
dds=subset(dds, select=(colData(dds)[, cond]%in%c(condA, condB)))
colData(dds)[, cond]=relevel(colData(dds)[, cond], ref=condB)
colData(dds)[, cond]=droplevels(colData(dds)[, cond])
print('==> h5ad2deseq2()')
str(dds)

if (rowmean>0) {
	dds=subset(dds, rowMeans(counts(dds, normalized=F))>=rowmean)
	dds=subset(dds, select=colSums(counts(dds, normalized=F))>0)
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

# Union
udds=DESeq(
	dds
	, fitType='glmGamPoi'
	, test='LRT'
	, useT=T
	, minmu=1e-6
	, minReplicatesForReplace=Inf
	, reduced=~1
	, parallel=F
	)
cat('==> DESeq(), union\n')
cat('    full: ~ group + cond\n')
cat('    reduced: ~ 1\n')
str(udds)
saveRDS(udds, file=sprintf('%s/%s_union_DESeq.rds', outdir, bname))
res=data.frame(results(udds))
res=res[order(res$padj), ]
write.txt(cbind(symbol=rownames(res), res), file=sprintf('%s/%s_union.txt.gz', outdir, bname))

# Consensus
cdds=DESeq(
	dds
	, fitType='glmGamPoi'
	, test='LRT'
	, useT=T
	, minmu=1e-6
	, minReplicatesForReplace=Inf
	, reduced=as.formula(sprintf('~%s', group))
	, parallel=F
	)
cat('==> DESeq(), consensus\n')
cat('    full: ~ group + cond\n')
cat('    reduced: ~ group\n')
str(cdds)
saveRDS(cdds, file=sprintf('%s/%s_consensus_DESeq.rds', outdir, bname))
res=data.frame(results(cdds))
res=res[order(res$padj), ]
write.txt(cbind(symbol=rownames(res), res), file=sprintf('%s/%s_consensus.txt.gz', outdir, bname))
