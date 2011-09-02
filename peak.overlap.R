# Justin Ashworth 2010
# jashworth@systemsbiology.org
# Institute for Systems Biology

peak.overlap =
	function(peaks,distance_cutoff=500,binstep=25,glength=2014239,matchdist=100,writefiles=T,...)
	# computes overlap of [ChIP] peaks
	# note default Halo genome length (used to generate empirical null distrubution)
{
	out=list()
	out[["table"]]=data.frame()
	peaksnames=names(peaks)
	l=length(peaksnames)
	out[["pairs"]]=matrix(nrow=l,ncol=l,dimnames=list(peaksnames,peaksnames))

	for(i in 1:l){
		aname=peaksnames[i]
		apos=peaks[[aname]]

		for(j in i:l){
			bname=peaksnames[j]
			pairkey=paste(aname,'.vs.',bname,sep='')
			print(pairkey)
			bpos=peaks[[bname]]
			bins=seq(0,distance_cutoff,binstep)
			pairing=pairwise.site.distance(apos,bpos)
			pairing=subset(pairing,dis<distance_cutoff)
			h=hist(pairing$dis,breaks=bins,plot=F)
			# null distribution: randomly select same numbers of sites in genome range
			nd=pairwise.site.distrib(length(apos),length(bpos),glength,cutoff=distance_cutoff,step=binstep,ntrials=50)

			if(writefiles){
				# write plot
				png(paste(pairkey,'.png',sep=''))
				plot(h,main=pairkey,xlab='site distance',...)
				mtext(paste(length(apos),'sites vs.',length(bpos),'sites'),line=-1)
				lines(nd)
				high=h$counts>(nd$mean+2*nd$sd)
				if(length(h$counts[high])>0) text(h$mids[high],h$counts[high]*0.8,labels=h$counts[high],srt=90)
				dev.off()
				# write bin table
				bintable=data.frame(
					binstart=bins[1:length(bins)-1],
					binend=bins[2:length(bins)],
					counts=h$counts,
					null=nd$mean,
					nullsd=nd$sd)
				write.table(bintable,paste(pairkey,'.table',sep=''),quote=F,sep='\t',row.names=F)
				# write pairings
				if(aname!=bname){
					write.table(
						round(subset(pairing,pairing$dis<matchdist)),
						paste(pairkey,'.overlap',sep=''),quote=F,sep='\t',row.names=F
					)
				}
			}

			# store 'significant' pairings
			if(aname!=bname){
				out[["pairs"]][[aname,bname]]=pairkey
				out[["pairs"]][[bname,aname]]=pairkey
				out[[pairkey]]=pairing
			}

			# tabulate close peaks for return value
			matchbins=h$mids<matchdist
			matches=sum(h$counts[matchbins])
			nullmatches=sum(round(nd$mean[matchbins]))
			out[["table"]]=
				rbind(out[["table"]],data.frame(chipA=aname,chipB=bname,matches=matches,nullmatches=nullmatches))
			# add symmetric table entry for convenience with later analysis that is not symmetrically aware
			# (e.g. xtabs/heatmap)
			if(aname!=bname)out[["table"]]=
				rbind(out[["table"]],data.frame(chipA=bname,chipB=aname,matches=matches,nullmatches=nullmatches))
		}
	}
	class(out)=c('peak.overlap',class(out))
	return(out)
}

library(gplots)
plot.peak.overlap =
	function(x,subtractnull=T,sep=c(),cellnotes=T,...)
	# plots heatmap of results from peak.overlap
{
	# list of pairs to 2D matrix using xtabs (note: asymmetric)
	# must specify factor levels explicitly to avoid annoying autosort
	x=x[["table"]]
	matches=x$matches
	if(subtractnull)matches=matches-x$nullmatches
	xt=xtabs(x$matches~factor(x$chipA,levels=unique(x$chipA))+factor(x$chipB,levels=unique(x$chipB)))
	if(nrow(xt)<2 | ncol(xt)<2){
		print('not enough data')
		return()
	}
	# must manually resize margins to fit axis labels
	op=par(oma=c(20,2,2,20))
#	efact=1.1
#	colors=rgb(efact^(32:1),efact^(32:1),efact^(32),maxColorValue=efact^(32))
	colors=rgb(32:1,32:1,32,maxColorValue=32)
	sepcol=rgb(0.6,1,0.6)
#	heatmap(xt,Rowv=NA,Colv="Rowv",cexRow=2,cexCol=2,scale='none',col=colors)
	cellnote=matrix(NA,nrow(xt),ncol(xt))
	# avoid labeling cells if too many
	if(cellnotes & length(cellnote)<400)cellnote=xt
	# Rowv, Colv arguments required to yield simple heatmap
	heatmap.2(xt,
		Rowv=NA,Colv="Rowv",key=F,
		scale='none',col=colors,colsep=sep,rowsep=sep,sepcolor=sepcol,
		trace='none',dendrogram='none',
		cellnote=cellnote,notecol='black',
		main='number of overlapping peaks',...)
	par(op)
}

pairwise.site.distance =
	function(asites,bsites)
{
	# set of all pairwise site combinations
	pairing=expand.grid(asites,bsites)
	names(pairing)=c('peak1','peak2')
	# distances
	pairing$dis=abs(pairing$peak2-pairing$peak1)
	invisible(pairing)
}

pairwise.site.distrib =
	# empirical distribution of distances between variable numbers of randomly chosen sites on a finite 1D axis
	# (e.g. ChIP peaks along genome coordinate)
	function(asites,bsites,l,cutoff=10000,step=100,ntrials=10)
{
	if(cutoff>l){cutoff=l}
	if(step>=cutoff){step=cutoff/10}
	bins=seq(0,cutoff,step)
	bincounts=matrix(rep(NA,cutoff/step*ntrials),ncol=ntrials)
	for(i in 1:ntrials){
		pairing=pairwise.site.distance( runif(asites,1,l), runif(bsites,1,l) )
		# below cutoff
		pairing=subset(pairing,dis<cutoff)
		h=hist(pairing$dis,breaks=bins,plot=F)
		bincounts[,i]=h$counts
	}

	value=data.frame(
		start=bins[1:length(bins)-1],
		end=bins[2:length(bins)],
		mean=rowMeans(bincounts),
		sd=sd(t(bincounts)))
	# R pseudo-object-orientation: assign additional custom class name to extend generic functions like 'plot' (below)
	class(value)=c('pairwise.site.distrib',class(value))
	invisible(value)
}

plot.pairwise.site.distrib =
	function(x,...)
{
	plot(0,0,xlim=c(min(x$start),max(x$end)),ylim=c(min(x$mean),max(x$mean)),type='n',...)
#	plot(x$h,...)
	lines(x,...)
}

lines.pairwise.site.distrib =
	function(x,...)
{
	mids=(x$start+x$end)/2
	lines(mids,x$mean,...)
	lines(mids,x$mean+x$sd,lty=2,...)
	lines(mids,x$mean-x$sd,lty=2,...)
}

filter_peaks =
	function(peaks,pval=0.05,trans_start_sites=NULL)
	# filters by p-value and by proximity to transcription start sites
{
	filtering=data.frame()
	for(id in names(peaks)){
		peakset=peaks[[id]]
		start=max(0,nrow(peakset))
		# column 1 contains the peak positions
		peakset=peakset[peakset[,3]<pval,1]
		passp=length(peakset)
		neartss=NA
		if( is.data.frame(trans_start_sites) ){
			peakset=peakset[sites_nearby(peakset,trans_start_sites$TSS)]
			neartss=length(peakset)
		}
		filtering=rbind(filtering,data.frame(peaks=start,passp=passp,neartss=neartss))
		peaks[[id]]=peakset
	}
	print(filtering)
	return(peaks)
}

sites_nearby =
	function(query,sites,cutoff=100)
{
#	return( query[ cutoff > sapply(query,function(x){min(abs(x-sites))}) ] )
	return ( cutoff > sapply(query,function(x){min(abs(x-sites))}) )
}

getseqFASTA =
	function(FASTA,key)
	# shortcut for extracting sequence(s) from Biostrings FASTA object (unwieldy list of pairs)
{
	select=lapply(FASTA,'[[','desc') %in% c(key)
	return( sapply(FASTA[select],'[[','seq') )
}

segmentsFASTA =
	function(startend,refseq,prefix=NULL)
{
	require(Biostrings)
	startend$seq=substring(refseq,startend$start,startend$end)
	# construct list of lists for writeFASTA
	l=lapply(	seq_len(nrow(startend)),
		function(x){ list(desc=paste(startend$start[x],'_',startend$end[x],sep=''),seq=startend$seq[x]) } )
	if( !is.null(prefix) ) writeFASTA(l,paste(prefix,'fa',sep='.'))
	else return(l)
}

medichi.fit.info =
	function(rdatafile)
{
	load(rdatafile)
	# median probe resolution
	coords=sort(as.numeric(fits$fits.fin$Chr[[1]]$data[,1]))
	d=diff(coords)
	m=median(d[d>1])
	# number of peaks
	npeaks=nrow(fits$fits.fin$Chr[[1]]$coeffs)
	return(list(resol=m,npeaks=npeaks))
}
