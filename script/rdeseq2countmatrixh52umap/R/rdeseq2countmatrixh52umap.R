# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(dplyr))

tmp=read.txt(label)
rownames(tmp)=tmp[, 1]
tmp[, 1]=NULL
print('==> label')
str(tmp)

x=Read10X_h5(f)
x=CreateSeuratObject(counts=x)

if (nrow(x)>ntop) {
	x=x %>% FindVariableFeatures(nfeatures=ntop) %>% ScaleData() %>% RunPCA() %>% RunUMAP(dims=1:30)
} else {
	x=x %>% ScaleData(features=rownames(x)) %>% RunPCA(features=rownames(x)) %>% RunUMAP(dims=1:30)
}

x=subset(x, cells=rownames(tmp))
x=AddMetaData(x, metadata=tmp)
print('==> x')
str(x)
write.txt(data.frame(barcode=rownames(x@meta.data), x@meta.data), file=sprintf('%s/%s_metadata.txt.gz', outdir, bname))
SeuratDisk::SaveH5Seurat(x, sprintf('%s/%s.h5seurat', outdir, bname))
seuratdimplotgroupby(x, reduct='umap', group=names(tmp)[1], label=F, outfile=sprintf('%s/%s.png', outdir, bname))
