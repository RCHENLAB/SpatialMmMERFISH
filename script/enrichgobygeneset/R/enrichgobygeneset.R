# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
qlist=read.txt(infile, header=F)[, 1, drop=T]
suppressPackageStartupMessages(library(clusterProfiler))
x=enrichGO(qlist, species, ont='all', qvalueCutoff=qthr, keyType='SYMBOL')
saveRDS(x, file=sprintf('%s/%s.rds', outdir, bname))
write.txt(x@result, file=sprintf('%s/%s.txt.gz', outdir, bname))
