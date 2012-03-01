ratiocolregex = 'Day'

ann_sig = function(ratios,sort=T) {
	upearly = read.delim('all.upearly')
	uplate  = read.delim('all.uplate')
	updark  = read.delim('all.updark')
	uplight = read.delim('all.uplight')
	ratios$early = ratios$geneid %in% upearly$geneid
	ratios$late  = ratios$geneid %in% uplate$geneid
	ratios$dark  = ratios$geneid %in% updark$geneid
	ratios$light = ratios$geneid %in% uplight$geneid
	if (sort) {
		ratiocols = grep(ratiocolregex,names(ratios))
		ord = order( apply(ratios[,ratiocols],1,sd), decreasing=T )
		ratios = ratios[ ord, ]
		ord = order( ratios$protein )
		ratios = ratios[ ord, ]
		hc = hclust( dist( ratios[,ratiocols] ) )
		ord = hc$order
		ratios = ratios[ ord, ]
		ratios = rbind(
			ratios[ ratios$early, ],
			ratios[ ratios$late, ],
			ratios[ ratios$dark, ],
			ratios[ ratios$light, ],
			ratios[ !ratios$early & !ratios$late & !ratios$dark & !ratios$light, ]
		)
	}
	return(ratios)
}

if (!exists('ratios')) {
	ratios = read.delim('ratios')
	genes = read.delim('chromatin.tsv')
	ratios = ratios[ ratios$geneid %in% genes$geneid, ]
	ratiosig = ann_sig(ratios)
	ratios = ratiosig
}

ratiocols = grep(ratiocolregex,names(ratios))
ratiocols = c(ratiocols,grep('early|late|dark|light',names(ratios)))
mat = as.matrix( ratios[,ratiocols] )
#rownames(mat) = ratios$name
#rownames(mat) = ratios$KOG.Description
#rownames(mat) = ratios$desc
rownames(mat) = sprintf('%-10s %-20s %s',ratios$geneid,ratios$tag,ratios$desc)
#mat = mat[1:18,]

library(gplots)
cols = colorpanel(64,'blue','black','yellow')

pdf('chromatin.pdf',width=20,height=15)
collabs = gsub('X','',colnames(mat))
collabs = gsub('\\.',' ',collabs)
collabs = gsub('Dk','Dark',collabs)
collabs = gsub('Lt','Light',collabs)
breaks = seq(-2,2,by=4/64)
op=par(oma=c(4,2,4,20),cex=2,cex.lab=2,cex.axis=2,cex.main=2,cex.sub=2)
cexRow=1.8
cexCol=cexRow
hm = heatmap.2(mat,trace='none',scale='none',Colv=NA,Rowv=NA,col=cols,lhei=c(1,5),lwid=c(4,3),density.info='none',cexRow=cexRow,cexCol=cexCol,labCol = collabs,breaks=breaks )
dev.off()
