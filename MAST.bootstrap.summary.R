logs = readLines('logs')

progress = lapply( logs,
	function(x){
		d = read.delim(x, skip=3, sep=' ', header=F, as.is=T)
		targets = do.call(rbind, strsplit( d[,1], '/'))
		background = do.call(rbind, strsplit( d[,9], '/'))
		pval = gsub('pval=','', d[,14])
		d = as.data.frame( cbind(targets, background, pval), stringsAsFactors=F )
		names(d) = c('targets.match','targets','all.match','all','pval')
		return(d)
	})
names(progress) = logs

start = do.call(rbind, lapply(progress, function(x){ x[1,] }))
final = do.call(rbind, lapply(progress, function(x){ x[nrow(x),] }))

names(start) = paste(names(start),'.start',sep='')
names(final) = paste(names(final),'.final',sep='')

start.final = cbind(start,final)
write.table(start.final,'logfile.table.tsv')

colkey = list(
	chip = '#0000FFBB',
	cor.chip = '#FF00FFBB'
)

cols = rep(colkey$chip,length(logs))
cols[ grepl('pos.neg',logs) ] = colkey$cor.chip

pdf('motif.match.progress.pdf',useDingbats=F)

plot(0,0,type='n', xlim=c(1,max(sapply(progress,nrow))), ylim=c(0,1), xlab='Rounds of fitting', ylab='fraction of genes w/ hits', xaxs='i', yaxs='i')
for(i in 1:length(progress)){
	prg = apply(progress[[i]],2,as.numeric)
	lines(1:nrow(prg),prg[,1]/prg[,2],col=cols[i],lwd=3)
}
legend('topleft',c('ChIP','correlated.ChIP'),lty=1,col=c(colkey$chip,colkey$cor.chip),lwd=3)
dev.off()

pdf('motif.pval.progress.pdf',useDingbats=F)

minp = as.numeric( min( unlist( sapply(progress,function(x){ as.numeric(x[,5]) })) ) )

plot(0,0,type='n', xlim=c(1,max(sapply(progress,nrow))), ylim=c(minp,1), xlab='Rounds of fitting', ylab='motif enrichment  p-value', xaxs='i', yaxs='i', log='y')
abline(h=0.05,lty=2,lwd=2)
for(i in 1:length(progress)){
	prg = apply(progress[[i]],2,as.numeric)
	lines(1:nrow(prg),prg[,5],col=cols[i],lwd=3)
}
legend('bottomleft',c('ChIP','correlated.ChIP'),lty=1,col=c(colkey$chip,colkey$cor.chip),lwd=3)
dev.off()
