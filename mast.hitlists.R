# for input files, this R script assumes that following has been done to format the hitlist as a tsv
#	echo "target	strand	start	end	score	pval" > mast.$lrp
#	[mast] -hit_list | sed -E 's/ +/	/g' | grep -v "#" >> hitlist
# e.g.:
#	target	strand	start	end	score	pval
#	VNG2604G	+1	249	263	-11.14	8.47e-04
#	VNG2604G	-1	365	379	-7.55	3.50e-04
#	VNG6404H	-1	320	334	-8.84	4.88e-04
#	VNG2146H	-1	232	246	-3.49	1.25e-04
#	VNG0966G	-1	189	203	4.25	1.34e-05
#	VNG6182H	+1	154	168	-3.03	1.11e-04
#	VNG0563G	+1	33	47	-6.26	2.55e-04
#	VNG1123G	+1	322	336	0.10	4.65e-05
#	VNG0146H	-1	157	171	-9.58	5.79e-04

# proper OO kind of approach would first translate all TF -> targetmotif lists into master gene -> TF -> targetmotif data structure

mast.hitlist.to.gene.centric.list =
	function(hitlist)
{
	genes = as.character(levels(factor(hitlist$target)))
	gl = lapply(genes,function(x){hitlist[hitlist$target==x,2:ncol(hitlist)]})
	names(gl) = genes
	class(gl)=c("mast.hitlist",class(gl))
	return(gl)
}

# load files, add columns to keep track of TFs after subsequent merging
load.mast.hitlists =
	function(pathsfile)
{
	paths=as.character(readLines(pathsfile))
	color=rainbow(length(paths),alpha=0.5,start=.25,end=.1)
	hitlists=lapply(1:length(paths),
		function(x){
			path=paths[x]
			hitlist=read.delim(path)
			hitlist$tf=gsub('mast.','',path)
			hitlist$color=color[x]
			return(hitlist)
		}
	)
	names(hitlists)=gsub('mast.','',paths)
	return(hitlists)
}

plot.mast.hitlist =
	function(x,xlim=c(-500,100),...)
{
	genes=names(x)
	ngenes=length(genes)
	plotindex=1
	op=par(mar=c(5,6,4,6))
	plot(NULL,xlim=xlim,ylim=c(0,ngenes+1),main='motif match positions',xlab='gene sequence position',ylab='genes',xaxs='i',yaxs='i',yaxt='n',...)
	abline(v=0)
	for(gene in genes){
		text(xlim[1],plotindex,gene,pos=4)
#		plot.gene.motifs(x[[gene]],xlim=xlim,plotindex=plotindex,name=gene,...)
		plot.gene.motifs(x[[gene]],xlim=xlim,plotindex=plotindex,...)
		plotindex=plotindex+1
	}
	# plot legend
	# identify and color the different TFs
	# this is ugly--is there a better way to do double-indexing in list of data.frames?
	colortab=unique(cbind(unlist(sapply(x,'[[','tf')),unlist(sapply(x,'[[','color'))))
	colortab=colortab[order(colortab[,1]),]
	mtext('TFs:',side=3,adj=0.99,line=-5,col='black',outer=T)
	for(i in 1:nrow(colortab)){
		mtext(colortab[i,1],side=3,adj=0.99,line=-1*i-5,col=colortab[i,2],outer=T)
	}
}

plot.gene.motifs =
	function(x,xlim=c(-500,100),plotindex=1,ntadj=0,name='',yscale=1,gthick=0.1,color='black')
	# 'x' is an object with x$strand, x$start, x$end, x$pval, x$color (potentially vectorized)
{
	# extension line
#	segments(x0=c(xlim[1]-200,xlim[2]),x1=c(xlim[1],xlim[2]+200),y0=rep(plotindex,4))
#	segments(x0=xlim[1]-200,x1=xlim[1]+ntadj,y0=rep(plotindex,2))
	# "gene"
	gthickscale=yscale*gthick
	# gene box
	rect(xleft=xlim[1]+ntadj,xright=xlim[2],ybottom=-1*gthickscale+plotindex,ytop=gthickscale+plotindex)
	topstrand = x$strand %in% c(1,'top','+','fwd')
	hmod = rep(1,length=length(x$strand))
	hmod[!topstrand]=-1
	# the motif box
	if(!is.null(x$color)) color=x$color
	rect(xleft=x$start+xlim[1]+ntadj,xright=x$end+xlim[1]+ntadj,ytop=(-1*log(x$pval,10)-2)/5*hmod+plotindex,ybottom=gthickscale*hmod+plotindex,col=color)
	text((xlim[2]-xlim[1])/2+xlim[1],plotindex,name)
}

plot.mast.hits.for.lrps.targets =
	function()
{
	hls=load.mast.hitlists('mast.hits')
	hlss=do.call(rbind,hls)
	gl=mast.hitlist.to.gene.centric.list(hlss)
	plot(gl,xlim=c(-600,100),ntadj=100)
	return(gl)
}

# could also cheat around data translation approach by indexing graphical positions for target genes and just running through TF->targetmotif lists
