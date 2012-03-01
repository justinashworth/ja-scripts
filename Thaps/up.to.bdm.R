early = read.delim('all.upearly')
late = read.delim('all.uplate')
dark = read.delim('all.updark')
light = read.delim('all.uplight')
all = read.delim('ratios.changesfromavg')
bdm = data.frame(
	gene=all$geneid,
	early=all$geneid %in% early$geneid,
	late=all$geneid %in% late$geneid,
	dark=all$geneid %in% dark$geneid,
	light=all$geneid %in% light$geneid)

bdm[,c(2:5)] = apply(bdm[,c(2:5)],2,as.integer)

tsv(bdm,'all.bdm')
