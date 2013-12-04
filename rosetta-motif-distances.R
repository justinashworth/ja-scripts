paths=readLines('dnas')
#paths=rep('2e1c.WT/design_mutations.dna',4)
bg=read.table('halo.bg.file',header=F)
bg=bg[,2]
cat('background probabilities:\n')
cat(bg)
cat('\n')
names=dirname(paths)
names=gsub('/design_mutations.dna','',names)
names=gsub('2e1c.P.OT3.FL11.','',names)
source('rosetta-motifs.R')
pdf=F
mats=lapply(1:length(paths),function(i){
	path=paths[i]
	name=names[i]
	cat(name,' ',path,'\n')
	if(pdf){pdf(paste(name,'.pdf',sep=''))}
	sl=seqlogoFromFasta(path,plot=pdf)
	if(pdf){dev.off()}
	return(sl) })
names(mats)=names

disfunc='ED'
#disfunc='KL'
#disfunc='ALLR'
dmats=ppm.dis.matrix(mats,disfunc=disfunc,bg=bg)

require(gplots)

oma=c(7,3,3,7)

pdf=T
cols=redblue(32)
if(pdf){pdf('rosetta.motif.distances.pdf');par(oma=oma)}else{dev.new();par(oma=oma)}
heatmap.2(dmats,symm=T,trace='none',scale='none',dendrogram='row',main='Rosetta motif distances',col=cols,symbreaks=F,symkey=F)
if(pdf){dev.off()}

# heatmap of just Halo TFs
if(pdf){pdf('rosetta.motif.distances.halo.pdf');par(oma=c(5,5,5,5))}else{dev.new();par(oma=oma)}
dmats.halo=dmats[grepl('Halo',rownames(dmats)),grepl('Halo',colnames(dmats))]
heatmap.2(dmats.halo,symm=T,trace='none',scale='none',dendrogram='row',main='Rosetta motif distances',col=cols,symbreaks=F,symkey=F)
if(pdf){dev.off()}
