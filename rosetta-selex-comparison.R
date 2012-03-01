source('rosetta-motifs.R')

paths.rosetta=readLines('dnas')
paths.selex=readLines('Yokoyama.SELEX.motifs')

bg=read.table('halo.bg.file',header=F)
bg=bg[,2]
cat('background probabilities:\n')
cat(bg)
cat('\n')

# load/format probability matrices (horizontal)

# rosetta position probability matrices (horizontal)
# filter rosetta names
names.rosetta=paths.rosetta
names.rosetta=gsub('/design_mutations.dna','',names.rosetta)
names.rosetta=gsub('2e1c.P.OT3.FL11.','',names.rosetta)
names.rosetta=basename(names.rosetta)
png=0
png=1
ppms.rosetta=lapply(1:length(paths.rosetta),function(i){
	path=paths.rosetta[i]
	name=names.rosetta[i]
	cat(name,' ',path,'\n')
	if(png){png(paste(name,'.png',sep=''),width=800,height=320)}
	ppm=seqlogoFromFasta(path,plot=png)
	if(png){dev.off()}
	return(ppm) })
names(ppms.rosetta)=names.rosetta

# read, convert selex matrices to compatible format
names.selex = basename(paths.selex)
names.selex = gsub('.motif.probs','',names.selex)
# tag with suffix
names.selex = paste(names.selex,'.SELEX',sep='')
ppms.selex=lapply(1:length(paths.selex),function(i){
	path=paths.selex[i]
	name=names.selex[i]
	cat(name,' ',path,'\n')
	probs=read.delim(path,comment.char='#')
	horiz=t(probs[,c(2:ncol(probs))])
	norm=t(t(horiz)*1/colSums(horiz))
	if(png){
		png(paste(name,'.png',sep=''),width=800,height=320)
		seqLogo(norm,ic.scale=F)
		dev.off()
	}
	return(norm) })
names(ppms.selex)=names.selex

# combine
all.ppms = c(ppms.rosetta,ppms.selex)

non.zero = 1e-2
# must avoid zeros for most distance metrics: add non.zero to all values and renormalize
all.ppms=lapply(all.ppms,function(x){x=x+non.zero;scale(x,center=F,colSums(x))})

# pairwise motif distance matrices
#disfunc='KL'
#disfunc='ALLR'
#disfunc='BLiC.inv'
disfunc='ED'
dmats=ppm.dis.matrix(all.ppms,disfunc=disfunc,bg=bg)

require(gplots)

oma=c(7,3,3,7)
op=par(oma=oma)

# various subsets (row/colnames)
rosetta = dmats[!grepl('SELEX',rownames(dmats)),!grepl('SELEX',colnames(dmats))]
selex_motifs = grep('SELEX',rownames(dmats),value=1)
selex_prefices = gsub('.SELEX','', selex_motifs )
selex_set = c(selex_motifs,selex_prefices)

dislab = disfunc
if (dislab == 'ED') dislab = 'Euclidian Distance'

main=paste('Motif distances\n','metric =',dislab)

#cols=redblue(64)
cols=colorpanel(64,'black','gray','white')

scale='none'
#scale='row'
oma=c(10,3,3,10)

if(png){png('rosetta.motif.distances.png',width=800,height=800);par(oma=oma)}else{dev.new();par(oma=oma)}
heatmap.2(rosetta,symm=T,trace='none',scale=scale,dendrogram='row',main=paste(main,'rosetta',sep='\n'),col=cols,symbreaks=F,symkey=F)
if(png){dev.off()}

# Rosetta: just Halo TFs
dmats.halo = rosetta[grepl('Halo',rownames(rosetta)),grepl('Halo',colnames(rosetta))]
if(png){png('rosetta.motif.distances.halo.png',width=800,height=800);par(oma=oma)}else{dev.new();par(oma=oma)}
heatmap.2(dmats.halo,symm=T,trace='none',scale=scale,dendrogram='row',main=paste(main,'rosetta(Halo)',sep='\n'),col=cols,symbreaks=F,symkey=F)
if(png){dev.off()}

# Rosetta: just SELEX proteins
rosetta_for_selex = rosetta[rownames(rosetta) %in% selex_prefices, colnames(rosetta) %in% selex_prefices]
if(png){png('motif.distances.pred.for.selex.png',width=800,height=800);par(oma=oma)}else{dev.new();par(oma=oma)}
heatmap.2(rosetta_for_selex,symm=T,dendrogram='row',trace='none',scale=scale,main=paste(main,'rosetta for SELEX',sep='\n'),col=cols,symbreaks=F,symkey=F,cellnote=round(rosetta_for_selex,2),notecol='white',notecex=1.3)
if(png){dev.off()}

# heatmap of SELEX all vs. all
#dmats.selex=dmats[rownames(dmats) %in% selex_set,colnames(dmats) %in% selex_set]
#dmats.selex = dmats.selex[!grepl('36',rownames(dmats.selex)),!grepl('36',colnames(dmats.selex))]
#dmats.selex = self.distance.to.NA(dmats.selex)
#if(png){png('motif.distances.selex.png',width=800,height=800);par(oma=oma)}else{dev.new();par(oma=oma)}
#heatmap.2(dmats.selex,symm=T,trace='none',scale='none',Rowv='none',dendrogram='col',main=main,col=cols,symbreaks=F,symkey=F)
#if(png){dev.off()}

# SELEX-SELEX (intra)
selex.intra = dmats[rownames(dmats) %in% selex_motifs, colnames(dmats) %in% selex_motifs]
if(png){png('motif.distances.selex.intra.all.png',width=800,height=800);par(oma=oma)}else{dev.new();par(oma=oma)}
heatmap.2(selex.intra,symm=T,dendrogram='row',trace='none',scale=scale,main=paste(main,'SELEX(intra)',sep='\n'),col=cols,symbreaks=F,symkey=F,cellnote=round(selex.intra,2),notecol='white',notecex=1.3)
if(png){dev.off()}

# SELEX-SELEX (intra) [excluding vars at pos 36]
selex.intra = selex.intra[!grepl('36',rownames(selex.intra)),!grepl('36',colnames(selex.intra))]
if(png){png('motif.distances.selex.intra.png',width=800,height=800);par(oma=oma)}else{dev.new();par(oma=oma)}
heatmap.2(selex.intra,symm=T,dendrogram='row',trace='none',scale=scale,main=paste(main,'SELEX(intra)',sep='\n'),col=cols,symbreaks=F,symkey=F,cellnote=round(selex.intra,2),notecol='white',notecex=1.3)
if(png){dev.off()}

# SELEX vs. Rosetta
selex.v.rosetta = dmats[rownames(dmats) %in% selex_prefices, colnames(dmats) %in% selex_motifs]
rownames(selex.v.rosetta)=paste(rownames(selex.v.rosetta),'.pred',sep='')
selex.v.rosetta = selex.v.rosetta[ order(rownames(selex.v.rosetta)), order(colnames(selex.v.rosetta)) ]
if(png){png('motif.distances.pred.vs.selex.all.png',width=800,height=800);par(oma=oma)}else{dev.new();par(oma=oma)}
heatmap.2(selex.v.rosetta,symm=T,Rowv=0,Colv='Rowv',dendrogram='none',trace='none',scale=scale,main=main,col=cols,symbreaks=F,symkey=F,cellnote=round(selex.v.rosetta,2),notecol='white',notecex=1.3)
if(png){dev.off()}

# heatmap of SELEX vs. Rosetta (excluding mutants at position 36)
selex.v.rosetta = selex.v.rosetta[!grepl('36',rownames(selex.v.rosetta)),!grepl('36',colnames(selex.v.rosetta))]
if(png){png('motif.distances.pred.vs.selex.png',width=800,height=800);par(oma=oma)}else{dev.new();par(oma=oma)}
heatmap.2(selex.v.rosetta,symm=T,Rowv=0,Colv='Rowv',dendrogram='none',trace='none',scale=scale,main=main,col=cols,symbreaks=F,symkey=F,cellnote=round(selex.v.rosetta,2),notecol='white',notecex=1.6,cexRow=2,cexCol=2,cex.main=2,key=FALSE)
#mtext('in vitro expt',2,adj=-0.2,cex=3)
if(png){dev.off()}

# clean version (for figure to be labeled manually)
rownames(selex.v.rosetta) = gsub( '.pred','',rownames(selex.v.rosetta) )
colnames(selex.v.rosetta) = gsub( '.SELEX','',colnames(selex.v.rosetta) )
if(png){png('motif.distances.pred.vs.selex.clean.png',width=800,height=800);par(oma=oma)}else{dev.new();par(oma=oma)}
heatmap.2(selex.v.rosetta,symm=T,Rowv=0,Colv='Rowv',dendrogram='none',trace='none',scale=scale,main=main,col=cols,symbreaks=F,symkey=F,cexRow=3,cexCol=3,cex.main=3,key=FALSE,cellnote=round(selex.v.rosetta,2),notecol='cyan',notecex=3)
#mtext('in vitro expt',2,adj=-0.2,cex=3)
if(png){dev.off()}

par=op
