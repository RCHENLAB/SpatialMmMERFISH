# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))

x=read.txt(f)
print('==> x')
str(x)

calproportion=function(tmp, save_plot=T, iter=0) {
	tmp=as.data.frame(table(tmp[, c(label, layer)]))
	tmp=dcast(tmp, formula(sprintf('%s~%s', label, layer)), value.var='Freq')
	tmp$ncell=rowSums(tmp[, -1])
	tmp$proportion=tmp[, group]/tmp$ncell
	tmp=tmp[order(-tmp[, group]), ]
	if (save_plot) {
		pdf(sprintf('%s/%s_prop_hist_iter%d.pdf', outdir, bname, iter), width=width, height=height)
		hist(tmp$proportion, xlab=sprintf('Proportion of cells displaced in %s', group), breaks=min(20, nrow(tmp)), main='', xlim=c(0, 1))
		dev.off()

		tmp[, label]=factor(tmp[, label], levels=unique(tmp[, label]))
		p=ggplot(data=tmp, aes_string(x=label, y='proportion')) +
		geom_bar(position='identity', stat='identity', show.legend=F, width=0.5, color='blue', fill='blue') +
		xlab(label) +
		ylab(sprintf('Proportion of cells displaced in %s', group)) +
		ggtitle('') +
		scale_y_continuous(breaks=scales::pretty_breaks(n=10), limits=c(0, 1), labels=scales::percent) +
		theme_bw() +
		theme(
			plot.background = element_blank()
			, panel.grid.major = element_blank()
			, panel.grid.minor = element_blank()
			, panel.border = element_blank()
			, plot.title = element_text(hjust=0.5)
			, axis.line = element_line(color='black')
			, axis.text.x = element_text(angle=90,vjust=0.5, hjust=1)
			, axis.text = element_text(color='black')
			)
		ggsave(p, file=sprintf('%s/%s_prop_%s_iter%d.pdf', outdir, bname, label, iter), width=width, height=height, units='in', useDingbats=F)
	}
	tmp
}
res=calproportion(x)
print('==> res')
str(res)

shufflebysection=function(dm, group, splitby) { # shuffle group within splitby
	do.call(
		rbind
		, by(
			dm
			, dm[[splitby]]
			, function(ds) {
				ds[[group]]=sample(ds[[group]])
				ds
			})
		)
}

iter=1
bgres=do.call(
	rbind
	, replicate(
		nreplicate
		, (
			function(){
				sx=shufflebysection(x, label, split)
				cat('.')
				sres=calproportion(sx, save_plot=F, iter)
				iter<<-iter+1
				sres
			}
			)()
		, simplify=F
		)
	)
write.txt(bgres, file=sprintf('%s/%s_bg.txt.gz', outdir, bname))

tmp=do.call(
	rbind
	, lapply(
		unique(res[[label]])
		, function(l) {
			obs=res[res[[label]]==l, 'proportion']
			bg=bgres[bgres[[label]]==l, 'proportion']
			pvalue=sum(bg>=obs)/length(bg)
			pdf(sprintf('%s/%s_%s_perm.pdf', outdir, bname, l), width=width, height=height)
			hist(bg, breaks=100, xlab='Proportion', main=sprintf('%s: obs=%s, pvalue=%G', l, scales::label_percent(accuracy=0.01)(obs), pvalue), xlim=c(min(c(obs, bg)), max(c(obs, bg))))
			abline(v=obs, col='red')
			dev.off()
			setNames(data.frame(l, pvalue), c(label, 'pvalue'))
		})
	)
res=merge(res, tmp, by=label)
res=res[order(res$pvalue, -res[, group]), ]
print('==> res')
print(res)
write.txt(res, file=sprintf('%s/%s.txt.gz', outdir, bname))
