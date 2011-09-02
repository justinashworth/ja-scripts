paths=readLines('dnas')
#paths=rep('2e1c.WT/design_mutations.dna',4)
bg=read.table('halo.bg.file',header=F)
bg=bg[,2]
cat('background probabilities:\n')
cat(bg)
cat('\n')
names=basename(paths)
names=gsub('/design_mutations.dna','',names)
names=gsub('2e1c.P.OT3.FL11.','',names)
source('rosetta-motifs.R')
png=1
mats=lapply(1:length(paths),function(i){
	path=paths[i]
	name=names[i]
	cat(name,' ',path,'\n')
	if(png){png(paste(name,'.png',sep=''),width=800,height=320)}
	sl=seqlogoFromFasta(path,plot=png)
	if(png){dev.off()}
	return(sl) })
names(mats)=names

disfunc='KL'
#disfunc='ALLR'
dmats=ppm.dis.matrix(mats,disfunc=disfunc,bg=bg)

require(gplots)

oma=c(7,3,3,7)

cols=redblue(32)
if(png){png('rosetta.motif.distances.png',width=800,height=800);par(oma=oma)}else{dev.new();par(oma=oma)}
heatmap.2(dmats,symm=T,trace='none',scale='none',dendrogram='row',main='Rosetta motif distances',col=cols,symbreaks=F,symkey=F)
if(png){dev.off()}

# heatmap of just Halo TFs
if(png){png('rosetta.motif.distances.halo.png',width=800,height=800);par(oma=c(8,8,8,8))}else{dev.new();par(oma=oma)}
dmats.halo=dmats[grepl('Halo',rownames(dmats)),grepl('Halo',colnames(dmats))]
heatmap.2(dmats.halo,symm=T,trace='none',scale='none',dendrogram='row',main='Rosetta motif distances',col=cols,symbreaks=F,symkey=F)
if(png){dev.off()}
