# vim: set noexpandtab tabstop=2:

suppressPackageStartupMessages(library(jlutils))
suppressPackageStartupMessages(library(pathfindR))

x=read.txt(infile)
print('==> x')
str(x)

result=run_pathfindR(
	x
	, gene_sets = gset
	, pin_name_path = pinpath
	, max_to_plot = ntop
	, convert2alias = symbolalias
	, n_processes = numthreads
	, output_dir = bname
	, silent_option = T
	)
print('==> result')
str(result)

write.txt(result, file=sprintf('%s.txt.gz', bname))
