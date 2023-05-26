# vim: set noexpandtab tabstop=2:

x=readRDS(infile)
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(ggplot2))
runontology=function(ontology) {
	xx=x %>% filter(ONTOLOGY==ontology)
	if (nrow(xx@result)>0) {
		p=cnetplot(xx, showCategory=min(ntop, nrow(xx@result)), cex_category=0.5, cex_label_category=0.5, cex_label_gene=0.5)
		ggsave(p, file=sprintf('%s/%s_%s_cnet.png', outdir, bname, ontology), width=width, height=height, units='in', dpi=500)

		p=heatplot(xx, showCategory=min(ntop, nrow(xx@result)))
		ggsave(p, file=sprintf('%s/%s_%s_heat.png', outdir, bname, ontology), width=width, height=height, units='in', dpi=500)

		xx=pairwise_termsim(xx)
		p=emapplot(xx, showCategory=min(ntop, nrow(xx@result)), cex_label_category=0.5, cex_category=0.5)
		ggsave(p, file=sprintf('%s/%s_%s_emap.png', outdir, bname, ontology), width=width, height=height, units='in', dpi=500)
	}
}
lapply(c('BP', 'MF', 'CC'), runontology)
