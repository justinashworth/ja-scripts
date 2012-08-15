multifile = FALSE
separate_files = TRUE
processfiles = multifile | separate_files

process.multifile =
	function(file)
{
	h=read.delim('hits')
	# TO DO HERE:
	# deal with -/+ prepended to motif index (hit direction) [previously done manually]

	# use motif index lookup table to add motif names to hits
	i=read.delim('../motif.index')
	names=i$motif[h$motif]
	h$motif_index = h$motif
	h$motif = names
	return(h)
}

process.separate_files =
	function(pathfile)
{
	files = readLines('hits.paths')
	hits = lapply(files, function(x){read.delim(x)} )
	names(hits) = gsub('.hits.tsv','',files)

	# separate file
	#hits file format:
	#sequence_name	motif	hit_start	hit_end	score	hit_p-value
	#Chr	+1	10145	10161	5.38	2.19e-05
	#Chr	+1	17850	17866	4.36	3.04e-05
	#Chr	-1	22800	22816	7.04	1.25e-05

	# the +/- in 'motif' column are hit direction
	# R strips the + signs upon reading
	# add separate 'dir' column to note direction
	for(x in names(hits)){
		hits[[x]]$dir='fwd'
		hits[[x]]$dir[grepl('-',hits[[x]]$motif)]='rvs'
	}

	# replace motif column with filename root
	for(x in names(hits)){
		hits[[x]]$motif=x
	}

	# collapse list of motif data frames to single data frame with rbind
	h = hits[[1]]
	for(i in 2:length(hits)){
		h = rbind(h,hits[[i]])
	}
	return(h)
}

# regions to annotate
if(!exists('regions')){
	upstream=200
	downstream=200
	coords=read.delim('TSS-TTS.txt')
	coords$TSS = as.numeric(coords$TSS)
	coords$cis.start = 0
	coords$cis.end = 0
	fwd = coords$Strand == 'For'
	coords$cis.start[fwd] = coords$TSS[fwd] - upstream
	coords$cis.end[fwd] = coords$TSS[fwd] + downstream
	rvs = coords$Strand == 'Rev'
	coords$cis.start[rvs] = coords$TSS[rvs] + upstream
	coords$cis.end[rvs] = coords$TSS[rvs] - downstream
	regions=coords[,names(coords) %in% c('Chr','Name','cis.start','cis.end')]
	names(regions)=c('seq','name','start','end')
	# filter out rows with NAs
	hasna = apply(regions, 1, function(x){ any(is.na(x)) } )
	regions = regions[ !hasna, ]
}

annotate.regions =
	# note: does not separate out sequence contigs
	function(posns,regns)
{
	return( sapply(posns,
		function(x){
			regns.contain = as.numeric(regns$start) <= x & as.numeric(regns$end) >= x
			if(all(! which(regns.contain))){regns.contain=NA}
			else{regns.contain=paste(regns$name[regns.contain],sep=',',collapse=',')}
			return(regns.contain)
		}
	))
}

gene.hit.table =
	# table: rows: regions, cols: motifs, value: region contains motif hit
	function(posns, regns)
{
	motifs = levels(factor(posns$motif))
	nmotifs = length(motifs)
	gene.hits = matrix(FALSE, nrow=nrow(regns), ncol=nmotifs, dimnames=list(regns$name,motifs))
	contigs = levels(factor(regns$seq))
	for(i in 1:nmotifs){
		motif = motifs[i]
		cat(motif,'\n')
		hits = posns[ posns$motif == motif, ]
		for(contig in contigs){
			cat(contig, '\n')
			incontig = as.character(hits$sequence_name) == contig
			cat(length(which(incontig)), 'in contig\n')
			hitpos = (hits$hit_start[incontig] + hits$hit_end[incontig])/2
			contighits =
				sapply(1:nrow(regns), function(x){
					regn = regns[x,]
					return( regn$seq==contig & any(as.numeric(hitpos) >= regn$start & as.numeric(hitpos) <= regn$end) )
				})
			cat(length(which(contighits)), 'hits\n')
			gene.hits[,i] = gene.hits[,i] | contighits
		}
	}
	return(gene.hits)
}

gene.hit.list =
	# list of lists: motifs->hits
	function(posns, regns)
{
	motifs = levels(factor(posns$motif))
	nmotifs = length(motifs)
	contigs = levels(factor(regns$seq))
	gene.hits = list()
	for(i in 1:nmotifs){
		motif = motifs[i]
		cat(motif,'\n')
		hits = posns[ posns$motif == motif, ]
		gene.hits[[motif]] = c()
		for(contig in contigs){
			cat(contig, '\n')
			incontig = as.character(hits$sequence_name) == contig
			cat(length(which(incontig)), 'in contig\n')
			hitpos = (hits$hit_start[incontig] + hits$hit_end[incontig])/2
			contighits =
				sapply(1:nrow(regns), function(x){
					regn = regns[x,]
					return( regn$seq==contig & any(as.numeric(hitpos) >= regn$start & as.numeric(hitpos) <= regn$end) )
				})
			genes = as.character(regns$name[contighits])
			gene.hits[[motif]] = unique( c(gene.hits[[motif]], genes) )
		}
	}
	return(gene.hits)
}

h = NULL
if (processfiles) {
	if (multifile) {
		path = 'hits'
		h = process.multifile(path)

	} else if (separate_files) {
		path = 'hits.paths'
		h = process.separate_files(path)
	}

	# rearrange columns
	# should yield:sequence_name	motif	dir	hit_start	hit_end	score	hit_p.value
	h = h[,c(1,2,7,3,4,5,6)]

#	h$regions = annotate.regions((h[,4]+h[,5])/2,regions)

#	gene.hits = gene.hit.table(h, regions)
#	write.table(gene.hits,'gene.hits.tab',quote=F,sep='\t')

	gene.hit.list = gene.hit.list(h, regions)

	# write file
	write.table(h,'hitsflat.tsv',quote=F,row.names=F,sep='\t')

} else {
	path = 'hitsflat.tsv'
	# this file should reflect everything done above
	h = read.delim(path)
}

# ANALYSIS/PLOTTING/OUTPUT
plots = FALSE
#plots = TRUE

if(plots){

png=FALSE
#png=TRUE

motifs = levels(factor(h$motif))
# filter down to Halo
motifs = motifs[grepl('Halo',motifs)]
# add back 2e1c.WT
#motifs = c(motifs,'2e1c.WT')
# filter hits
h = subset(h,motif %in% motifs)
require(gdata)
h=drop.levels(h)

if(png){png('genome.wide.hit.distrib.png',w=640,h=640)}
par(mar=c(5,10,5,2))
boxplot(log(h$hit_p.value,10)~h$motif,xlab='MAST motif hits log10(p-value)',las=2,horizontal=T,main='homology-based motif matches in whole H. salinarum genome')
if(png){dev.off()}

p.val = 1e-5
#cat('*** APPLYING MAST P-VALUE CUTOFF OF',p.val,'***\n')
#h=subset(h,hit_p.value < p.val)

counts=table(h$in.region,h$motif)
counts2=counts[,order(counts[2,],decreasing=T)]

if(png){png('motif.hits.in.regions.png',w=640,h=640)}
par(mar=c(10,10,10,6))
barplot(counts2,las=2,legend=c('not in region','in region'),col=c('darkblue','red'),ylab='number of motif hits',main='rosetta motif hits in Halo regions')
if(png){dev.off()}

p.enrich=counts2[2,]/(counts[1,]+counts2[2,])
p.enrich=sort(p.enrich,decreas=T)

if(png){png('motif.hits.in.regions.fraction.png',w=640,h=640)}
par(mar=c(10,10,10,6))
barplot(p.enrich,las=2,ylab='fraction of motif hits',main='fraction of rosetta motif hits in Halo regions\n(fraction of all genome hits per TF)')
if(png){dev.off()}

if(png){png('motif.hits.in.regions.p.values.png',w=640,h=640)}
par(mar=c(10,10,10,6))
cols=c('darkblue','red')
boxplot(log(h$hit_p.value,10) ~ h$in.region * h$motif,las=2,col=cols,outcol=cols,outpch=21,legend=c('not in region','in region'), horizontal=T,main='p-values for motifs in regions vs. not in regions')
if(png){dev.off()}
} # plots
