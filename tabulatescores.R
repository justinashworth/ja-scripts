tabulate.scores =
	function(scorefile, wtfile, suffix='')
{
	scores = read.delim(scorefile)
	wt = read.delim(wtfile)
	mw = merge(wt,scores,by.x=c('pos','wt'),by.y=c('pos','mut'))
	names(mw)[3] = 'bound_wt'
	m = merge(mw,scores,by.x='pos',by.y='pos')
	m$ddG = round(m$bound - m$bound_wt, 5)
	tsv(m, paste('ddGs',suffix,sep='.'))
	ndestab = sapply( tapply( m$ddG, m$pos, function(x){ x > 3 } ), function(x) { length(which(x)) } )
	write.table(ndestab, paste('ndestab',suffix,sep='.'), quote=F)
	return(m)
}
