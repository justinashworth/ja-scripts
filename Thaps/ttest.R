# Justin Ashworth
# Institute for Systems Biology
# 2011

ratiocolregex = 'Day'

do.ttest =
	function(ratiosfile='probe.ratios.mean.tsv')
{
	d = read.delim(ratiosfile)
	ratiocols = grep(ratiocolregex,names(d),value=T)
	ratios = as.matrix( d[ , ratiocols ] )
	colnames(ratios) = ratiocols
	day1 = grep('Day1',colnames(ratios))
	day4 = grep('Day4',colnames(ratios))
	d$pvs  =  apply( ratios, 1, function(x) { res = t.test( x[day1], x[day4] ); res$p.value } )
	d$p.adj = p.adjust(d$pvs)
	d$change = apply( ratios, 1, function(x) { mean( as.numeric(x[day4]) ) } )
	d$sd = apply( ratios, 1, function(x) { sd( as.numeric(x[day4]) ) } )
	d = d[ order(d$pvs), ]
	d = d[ , c( which( !grepl('Day',names(d)) ), grep('Day',names(d)) ) ]
	tsv( d, paste( gsub('.tsv','',ratiosfile), 'ttest.tsv', sep='.' ) )
	return(d)
}

plot.ttest =
	function(d)
{
	library(gplots)

	dd=d[ 1:80, ]
	ratiocols = grep(ratiocolregex,names(dd))
	sig.mat = as.matrix( dd[ , ratiocols ] )
	rownames(sig.mat) = dd$geneid

	png('ttest.day1.vs.day4.png',width=1000,height=1000)
	op=par(oma=c(12,0,0,16))
	heatmap.2( sig.mat, trace='none', scale='none', col=colorpanel(32,'blue','black','yellow'), dendrogram='none', labRow='', Rowv='none', Colv='none', main = 'Gwenn\'s chemostats, CO2 increasing\nttest, Day 1 vs. Day 4' )
	dev.off()

	# known genes
	ddd=dd[ !is.na(dd$desc) & !dd$desc=='' & !grepl('hypothetical',dd$desc), ]
	ratiocols = grep(ratiocolregex,names(ddd))
	sig.mat = as.matrix( ddd[ , ratiocols ] )
	row.names(sig.mat) = ddd$desc

	png('ttest.day1.vs.day4.known.png',width=1000,height=1000)
	op=par(oma=c(12,0,0,16))
	heatmap.2( sig.mat, trace='none', col=colorpanel(32,'blue','black','yellow'), dendrogram='none', Colv='none', main = 'Gwenn\'s chemostats, CO2 increasing\nt-test, Day 1 vs. Day 4' )
	dev.off()

#	d=read.delim('probe.ratios.mean.ttest.tsv')
	day1day4 = c( rep('Day1',4), rep('Day4',4) )
	day1 = grep( 'Day1', colnames(d) )
	day4 = grep( 'Day4', colnames(d) )
	png('ttest.day1.vs.day4.top.png',width=1000,height=1000)
	par(mfrow=c(3,3))
	par(cex=1.5)
	par(cex.main=0.75)
	apply( d[1:9,], 1, function(x) { stripchart( as.numeric(c(x[day1],x[day4]))~day1day4, main=paste('p.adj=',x[colnames(d)=='p.adj'],sep=''), ,vertical=T, xlab=x[1], ylab='log ratio',pch=21 ) } )
	dev.off()
}

#d = do.ttest()
#plot.ttest(d)
