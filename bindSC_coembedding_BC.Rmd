library(Seurat)
library(bindSC)
library(umap)
library(Matrix)
library(rdist)
library(ggplot2)
library(ggpubr)
library(e1071)
library(progress)
library(irlba)

## load scRNA-seq reference
RNA<-SeuratDisk::LoadH5Seurat("/storage/singlecell/jinli/wkfl/metaanalysis/mouse/Josh/integration/annotation/resolution2/subclustering/BC_rmdoublets/anno/scrnaseurath5ad2h5seurat/mm_BC.h5seurat", assay = 'RNA')

## load MERFISH data
ST<-SeuratDisk::LoadH5Seurat("MERFISH_cpmsintercept_BC.h5seurat", assay = 'RNA')

## run integration
geneOverlap<-intersect(rownames(ST), rownames(RNA))

X <-ST[["RNA"]][geneOverlap,]
Z <- RNA[["RNA"]][geneOverlap,]
Y <- RNA@reductions$scVI@cell.embeddings

x.clst<-ST@meta.data$celltype_scVI
y.clst <- RNA@meta.data$celltype

res <- BiCCA( X = X[geneOverlap,] ,
            Y = t(Y),
            Z0 =Z[geneOverlap,],
            X.clst = x.clst,
            Y.clst = y.clst,
            K = 15,
            alpha = 0.2,
            lambda = 0.2,
            temp.path  ="out_k15_al02_noblock",
            num.iteration = 50,
            tolerance = 0.01,
            save = TRUE,
            parameter.optimize = FALSE,
            calc.score=FALSE,
            block.size = 0)   

saveRDS(res, file = "out_k15_al02_noblock.rds")

x_lst<-seq(1,length(x.clst),1)
y_lst<-seq(1,length(y.clst),1)

x.clst<-as.character(x.clst)
y.clst<-as.character(y.clst)


## generate co-embeding UMAP
bindsc_umap <-  umap(rbind(res$u, res$r))

umap_plt <- bindsc_umap
umap_plt  <- data.frame("UMAP1"=umap_plt$layout[,1],
                      "UMAP2"=umap_plt$layout[,2],
                      "celltype" = as.character(c(x.clst[x_lst], y.clst[y_lst])),
                      "data" = c(rep("MERFISH",length(x.clst[x_lst])),
                                 rep("scRNA",length(y.clst[y_lst]))))
write.table(umap_plt, file="umap_plt_BC_k15_al02_noblock.txt")


## plot co-embeding UMAP
xlim <- c(min(umap_plt$UMAP1), max(umap_plt$UMAP1))
ylim <-c(-11,11)

p21 <- UMAP_plot(meta = umap_plt[umap_plt$data=="MERFISH",],mylabel=paletteDiscrete(umap_plt$celltype),
               color = "celltype", xlim = xlim, ylim = ylim,alpha = 1) + ggtitle("MERFISH")
p22 <- UMAP_plot(meta = umap_plt[umap_plt$data=="scRNA",], mylabel=paletteDiscrete(umap_plt$celltype),
               color = "celltype", xlim = xlim, ylim = ylim,alpha = 1)  + ggtitle("scRNA")

p3 <- ggarrange(p21,p22)

tiff("coembedding_BC_k15_al02_noblock.tiff", width=13, height =6, res =300, units = "in", compression = "lzw")
print(p3)
dev.off()


## label transfer to predict subtype annotation
xx=label_transfer(dt1=res$r, X.clst=y.clst, dt2=res$u)
write.table(data.frame(barcode=rownames(xx), xx), file="merfish_anno_BC_k15_al02_noblock.txt")


