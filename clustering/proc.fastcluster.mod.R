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
hc$method = 'fastcluster.average'
hc$call = 'fastcluster.bootstrap.cpp'
hc$dist.method = 'Pearson distance (correlation)'

# pick some reasonable cut heights based on the distribution
meanh = mean(hc$height)
sdh = sd(hc$height)
cat('mean height is:',meanh,'(sd',sdh,')\n')
# reasonable range of cutheights for 'average' method
hs = seq(meanh+sdh, max(hc$height)-2*sdh, sdh/5)

# reasonable range of cutheights for Ward method
loghs = seq(-0.5,1,0.1)
hs = 10^loghs

library(dendextendRcpp)

if(!exists('cuts')){
cuts = lapply(hs, function(h){
	cat('Cutting at height',h,':')
	cls = Rcpp_cut_lower(as.dendrogram(hc), h=h)
	cat(length(cls),'clusters\n')
	cls
})
names(cuts) = sapply(cuts,function(x){length(x)})
save(cuts,file='cuts.RData')
}

ann = read.delim('~/diatoms/genemodels/tps.all.models.ja.tsv',as.is=T)
desc = ann$desc
names(desc) = ann$id

# plot clusters containing specific genes
qs = c('233','262258','39799','35133','ThpsCp008','ThpsCp009','bd2088','268546','34543')
for(q in qs){
	cat('plotting clusters for',q,'\n')
	hasit = sapply(cuts, function(cls){ which( sapply(cls, function(x){q %in% labels(x)})) })
	names(hasit) = names(cuts)

	# plot subdendrogram
	for(i in 1:length(cuts)){
		if(length(hasit)==0) next
		if(length(hasit[[i]])==0) next
		if(is.na(hasit[[i]] | is.null(hasit[[i]]))) next
		cl = cuts[[i]][[hasit[[i]]]]
		labels(cl) = paste( labels(cl), desc[ labels(cl) ] )

		pdf(sprintf('cluster.%05i.%s.pdf',as.integer(names(cuts)[i]),q),useDingbats=F,height=1.5+length(unlist(cl))/6,width=10)
		par(mar=c(4,2,2,16))
		plot( cl, horiz=T, main=sprintf('cluster %i; cutheight=%g k=%s',hasit[[i]],hs[i],names(cuts)[i]))
		dev.off()
	}
}

clids = lapply(cuts, function(cut){ lapply(cut, function(x){labels(x)}) })
names(clids) = names(cuts)

for(k in names(clids)){
	pdf(sprintf('sizedist.%s.pdf',k),useDingbats=F)
	plot(density(sapply(clids[[k]],length)),main=sprintf('cluster size distribution: k=%s',k))
	dev.off()
}

write=F
if(write){
# write out clusters
choice = "651"
cat('cluster ids for',choice,'clusters\n')

if(!exists('rr')){
# load the corresponding NON-IMPUTED ratios
rr = read.table('ratios.lessthanhalfNAs.tab')
rr = as.matrix(rr)
}

cat('cluster stats...\n')
source('geneexp.tps.co2.lib.R')
cluster_stats(rr[unlist(clids[[choice]]),],clids[[choice]],spearman_means=T)
}
