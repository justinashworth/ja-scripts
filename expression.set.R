expression.set =
	function(ratios)
	# ratios is a 2D matrix with genes in rows and conditions in columns
{
	result=list()
	result$ratios=ratios
	# cor works on columns of matrices: use t() to compute table of all row correlations
	result$cor=cor(t(ratios))
	result$cor=result$cor[order(rownames(result$cor)),order(colnames(result$cor))]
	# 'OO magic' links result 'class' to general plot method below
	class(result)=c('expression.set',class(result))
	return(result)
}

plot.expression.set =
	function(x,new.dev=T,...)
{
	plot.expression.set.cor(x,new.dev,...)
	plot.expression.set.hc(x,new.dev,...)
}

plot.expression.set.cor =
	function(x,new.dev=T,...)
	# plot pairwise correlations
{
	require(gplots)
	if(new.dev) dev.new()
	op=par(oma=c(4,1,1,4))
	bins=seq(-1,1,0.1)
	colors=colorpanel(20,'red','white','blue')
	heatmap.2(x$cor,trace="none",cellnote=round(x$cor,2),scale="none",Colv="Rowv",col=colors,breaks=bins,notecol='black',key=F,lhei=c(0.1,1),...)
	par(op)
}

plot.expression.set.hc =
	function(x,new.dev=T,binrange=NULL,...)
	# plot hc heatmap of genes vs. conditions
{
	require(gplots)
	if(new.dev) dev.new(width=16)
	# filter out conditions with low variation
#	sigratios=x$ratios[,sd(x$ratios)>1]
	sigratios=x$ratios
	# compute distances
	hr=hclust(as.dist(1-cor(t(sigratios))))
	hc=hclust(as.dist(1-cor(sigratios)))
	op=par(oma=c(6,1,1,12))
	# make data bins for coloring
	nbins=20
	if(is.null(binrange)){
		maxabs=max(max(sigratios),abs(min(sigratios)))
		binrange=c(-maxabs,maxabs)
	}
	bins=seq(binrange[1],binrange[2],(binrange[2]-binrange[1])/nbins)
	colors=colorpanel(nbins,'blue','black','yellow')
#	heatmap.2(as.matrix(sigratios), Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=colors, breaks=bins, scale="none", trace="none",key=T,cexCol=1,lwid=c(0.1,1),...)
	heatmap.2(as.matrix(sigratios), Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=colors, breaks=bins, scale="none", trace="none",key=F,cexCol=2,cexRow=2,lwid=c(0.2,1),lhei=c(0.2,1),...)
#	rownames=hr$labels
#	rownames=lapply(hr$labels,function(x){do.call('expression',list(sprintf('Delta~%s',x)))})
#	heatmap.2(as.matrix(sigratios), Rowv=as.dendrogram(hr), Colv=as.dendrogram(hc), col=colors, breaks=bins, scale="none", trace="none",key=F,cexCol=1.5,cexRow=2,lwid=c(0.2,1),lhei=c(0.2,1),labRow=rownames,...)
	par(op)
}

write.expression.set =
	function(x,prefix='expression.set',...)
{
	# plots
	require(Cairo)
	Cairo(sprintf('%s.cor.png',prefix),type='png',width=720,height=600)
	plot.expression.set.cor(x,new.dev=F,...)
	dev.off()
	Cairo(sprintf('%s.hc.png',prefix),type='png',width=1000,height=320)
	plot.expression.set.hc(x,new.dev=F,...)
	dev.off()
}

halo.expression.set =
	function(genes)
{
	ratios=NULL
	if(file.exists('halo.ratios.tsv')) ratios=read.delim('halo.ratios.tsv')
	else {
		library(cMonkey)
		data(halo)
		ratios=halo$ratios
	}
	set.ratios=ratios[rownames(ratios) %in% genes$id,]
	# replace id with gene name
	names=as.character(sapply(rownames(set.ratios),function(x){genes$name[genes$id==x]}))
	rownames(set.ratios)=names
	return(expression.set(set.ratios))
}

demo.trh =
	function()
	# trh example
{
	trhnames=NULL
	if(file.exists('trh.codes')) trhnames=read.delim('trh.codes')
	else {
		trhnames=data.frame(
			id=c('VNG1404G','VNG1285G','VNG1816G','VNG2094G','VNG1922G','VNG1351G','VNG1123G','VNG1377G','VNG1179C','VNG1237C'),
			name=c('trh1','trh2','trh3','trh4','trh5','trh6','trh7','asnC','VNG1179C','VNG1237C')
		)
	}
	h=halo.expression.set(trhnames)
	plot(h)
	write.expression.set(h,prefix='trh')
	invisible(h)
}
