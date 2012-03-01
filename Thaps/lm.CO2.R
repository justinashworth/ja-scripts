# Justin Ashworth
# Institute for Systems Biology
# 2011

# CO2 levels for Gwenn's samples
# note assumed column order
day1 = 1:4
day2 = 5:8
day3 = 9:12
day4 = 13:16

co2 = rep(0,16)
#co2[day1]=600
#co2[day2]=900
#co2[day3]=1150
#co2[day4]=1135

co2 = c(
	316,
	310,
	330,
	322,

	372,
	362,
	399,
	393,

	464,
	458,
	454,
	451,

	553,
	455,
	433,
	441
)

ratiocolregex = 'Day'

probe.lm.cont =
	function( ratios, indep.var, plot=FALSE, ... )
	# to be called on a single probe/gene row, such as with apply()
{
	ratioframe = data.frame( indep = indep.var, dep = ratios )
	fit = summary( lm( dep ~ indep, ratioframe ) )
	if ( plot ) {
		plot(ratioframe,...)
		abline(fit)
	}
	result = c( fit$coefficients[2], pf( fit$fstatistic[1], fit$fstatistic[2], fit$fstatistic[3], lower.tail = F ) )
	names(result) = c('slope','pval')
	return(result)
}

lm.co2 =
	function(ratiosfile = 'probe.ratios.mean.tsv',unlog=FALSE,logbase=2)
{
	d = read.delim(ratiosfile)
	ratiocols = grep( ratiocolregex, names(d), value=T )
	if (unlog) d[,ratiocols] = logbase^d[,ratiocols]
	ratios = as.matrix( d[,ratiocols] )
	colnames(ratios) = ratiocols
	d = cbind( d, t( apply( ratios, 1, function(x) { probe.lm.cont(x,co2) } ) ) )
	d$p.adj = p.adjust(d$pval,'BH')
	d = d[ order(d$pval), ]
	d = d[ , c( which( !grepl(ratiocolregex,names(d)) ), grep(ratiocolregex,names(d)) ) ]
	tsv( d, paste( gsub('.tsv','',ratiosfile), 'lmCO3.tsv', sep='.' ) )
	return(d)
}

plot.lm.co2 =
	function(d)
{
	png('most.linear.with.co2.png',width=1000,height=1000)
	par(mfrow=c(3,3))
	par(cex=1.5)
	par(cex.main=0.75)
	apply( d[ 1:9, c( which(names(d) == 'geneid'),which(names(d) == 'desc'), which(names(d) == 'p.adj'), grep(ratiocolregex,names(d) ) ) ], 1,
#		function(x) { probe.lm.cont( as.numeric(x[4:length(x)]), co2, plot=TRUE, xlab='[CO2] (uatm)', ylab='log2 change from Day 1',main=paste(x[1],x[2],paste('p.adj=',x[3],sep=''),sep='\n') ) } )
		function(x) { probe.lm.cont( as.numeric(x[4:length(x)]), co2, plot=TRUE, xlab='[CO2] (uatm)', ylab='fold change',main=paste(x[1],x[2],paste('p(BH)=',x[3],sep=''),sep='\n') ) } )
	dev.off()

	png('exp.vs.lmCO2.png',width=800,height=800)
	par(oma=c(2,1,1,8))
	library(gplots)
	top = as.matrix( d[1:80,grep('Day',names(d))] )
	rownames(top) = d$geneid[1:80]
	heatmap.2( top, col=colorpanel(32,'blue','black','yellow'), Colv='none', Rowv='none', scale='none', dendrogram='none', trace='none', main='Gwenn\'s chemostats, CO2 increasing\ngenes most linear vs. CO2' )
	dev.off()
}

#d = lm.co2()
#plot.lm.co2(d)
