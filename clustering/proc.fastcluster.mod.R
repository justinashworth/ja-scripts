# script converts output from c++ fastcluster extension into R hclust object and dendrograms
# can probably re-implement the height-based tree cutting in c++ for bootstrapping/empirical distributions,
# and/or adapt code from source of dendextendRcpp's Rcpp_cut_lower function

# make a little dummy hclust object in R to fill with real data
hc = hclust(as.dist(matrix(runif(100),nrow=10)))
hc$merge = as.numeric( unlist( strsplit(readLines('hc.merges')[1],' ') ))
dim(hc$merge) = c( length(hc$merge)/2, 2 )
hc$height = as.numeric( unlist( strsplit(readLines('hc.heights')[1],' ') ))
hc$order = as.numeric( unlist( strsplit(readLines('hc.orders')[1],' ') ))
hc$labels = unlist( strsplit(readLines('hc.labels')[1],' ') )
hc$labels = gsub('\"','',hc$labels)
hc$method = 'fastcluster.average'
hc$call = 'fastcluster.bootstrap.cpp'
hc$dist.method = 'Pearson distance (correlation)'

cat('reading cached results (if any)\n')
if(file.exists('clusters.RData')) load('clusters.RData')
if(file.exists('hpk.RData')) load('hpk.RData')

library(dendextendRcpp)
if(!exists('hpk')){
	cat('Heights per k...\n')
	hpk = dendextendRcpp_heights_per_k.dendrogram(as.dendrogram(hc))
	save(hpk,file='hpk.RData')
}

ks = seq(100,2000,100)
hs = hpk[ as.character(ks) ]
#names(hs) = ks
#cat('heights:',hs,'\n')

if(!exists('cuts')){
cuts = list()
for(h in hs){
	cat('Cutting at height',h,':')
	cls = Rcpp_cut_lower(as.dendrogram(hc), h=h)
	l = length(cls)
	cat(l,'clusters\n')
	cuts[[as.character(l)]] = cls
	if(l==1) break
}
clids = lapply(cuts, function(cut){ lapply(cut, function(x){labels(x)}) })
names(clids) = names(cuts)
save.image(file='clusters.RData')
}

ann = read.delim('~/diatoms/genemodels/tps.all.models.ja.tsv',as.is=T)
desc = ann$desc
names(desc) = ann$id

if(!exists('rr')){
# load the corresponding NON-IMPUTED ratios
#rr = read.table('ratios.lessthanhalfNAs.tab')
cat('reading ratios\n')
rr = read.table('ratios')
rr = as.matrix(rr)
}

clplots=T
maxplot=300
if(clplots){
# plot clusters containing specific genes
qs = c('233','262258','39799','35133','268546','34543','33018')
#qs = as.character( c("233") )
for(q in qs){
	cat('plotting clusters for',q,'\n')
	hasit = sapply(cuts, function(cls){ which( sapply(cls, function(x){q %in% labels(x)})) })
	names(hasit) = names(cuts)

	# plot subdendrogram
	for(i in 1:length(cuts)){
		fname = sprintf('cluster.%05i.%s.pdf',as.integer(names(cuts)[i]),q)
#		if(file.exists(fname)){ cat('plot',fname,'exists, skipping\n'); next }
		if(length(hasit)==0) next
		if(length(hasit[[i]])==0) next
		if(is.na(hasit[[i]] | is.null(hasit[[i]]))) next
		cl = cuts[[i]][[hasit[[i]]]]
		sz = length(labels(cl))
		if(sz>maxplot) next
		ids = labels(cl)
		labels(cl) = paste( labels(cl), desc[ labels(cl) ] )

		baseh = 6
		pdf(fname,useDingbats=F,height=baseh+max(baseh/2,sz/4),width=baseh*1.5)
#		dev.new()
		#layout( matrix( c(1,2), nrow=2 ), widths=c(1,1), heights=c(1,sz/6))
		layout( matrix( c(1,2), nrow=2, ncol=1 ), heights=c(baseh,max(baseh/2,sz/4)), T )
		matplot( t(rr[ids,]), type='l', lwd=2, lty=1, col=rgb(0,0,0,0.5), main=fname )
		mtext(sprintf('%i genes',sz),side=3,line=-2)
		par(mar=c(4,2,2,10))
		plot(cl, horiz=T, main=sprintf('cluster %i; cutheight=%g k=%s',hasit[[i]],hs[i],names(cuts)[i]))
		dev.off()

	}
}

for(k in names(clids)){
	fname=sprintf('sizedist.%s.pdf',k)
#	if(file.exists(fname)){ cat('plot',fname,'exists, skipping\n'); next }
	pdf(fname,useDingbats=F)
	plot(density(sapply(clids[[k]],length)),main=sprintf('cluster size distribution: k=%s',k))
	dev.off()
}
}


choice = "500"

# plot all cluster dendrograms
cat('dendrograms\n')
maxdend=400
dir.create('dendro')
for(i in 1:length(cuts[[choice]])){
	cat('.')
	fname = sprintf('dendro/dendro.%i.png',i)
	cl = cuts[[choice]][[i]]
	ids = labels(cl)
	sz = length(ids)
	if(sz > maxdend) next
	labels(cl) = paste( labels(cl), desc[ labels(cl) ] )

	baseh = 480
	png(fname,height=max(baseh*0.5,sz*15),width=baseh)
	par(mar=c(2,1,2,20))
	plot(cl, horiz=T, main=sprintf('cluster %i',i))
	dev.off()
}
cat('\n')

write=T
if(write){
# write out clusters
cat('cluster stats for',choice,'clusters\n')
source('geneexp.tps.co2.lib.R')
cluster_stats(rr[unlist(clids[[choice]]),],clids[[choice]],spearman_means=T,co2rel=T)
}
