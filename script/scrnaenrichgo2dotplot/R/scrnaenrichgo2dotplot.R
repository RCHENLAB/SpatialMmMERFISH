# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(ggplot2))

x=readRDS(infile)
str(x)

runontology=function(ontology) {
	xx=x %>% filter(ONTOLOGY==ontology)
	if (nrow(xx@result)>0) {
		p=barplot(xx, showCategory=min(ntop, nrow(xx@result)), font.size=fontsize)
		outfile=sprintf('%s/%s_%s_barplot.pdf', outdir, bname, ontology)
		ggsave(p, file=outfile, width=width, height=height, units='in', useDingbats=F)

		p=dotplot(xx, showCategory=min(ntop, nrow(xx@result)), font.size=fontsize)
		outfile=sprintf('%s/%s_%s_dotplot.pdf', outdir, bname, ontology)
		ggsave(p, file=outfile, width=width, height=height, units='in', useDingbats=F)
	}
}
lapply(c('BP', 'MF', 'CC'), runontology)
