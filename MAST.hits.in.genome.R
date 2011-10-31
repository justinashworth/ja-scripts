multifile = 0
separate_files = 1
processfiles = multifile || separate_files

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

# promoter regions
if(!exists('regions')){
	upstream=500
	downstream=100
	coords=read.delim('halo.gene.coords.tsv')
	coords$Start = as.numeric(coords$Start)
	coords$cis.start = 0
	coords$cis.end = 0
	fwd = coords$Orientation == 'For'
	coords$cis.start[fwd] = coords$Start[fwd] - upstream
	coords$cis.end[fwd] = coords$Start[fwd] + downstream
	# in halo gene coords, 'Start' (generally) refers to actual start (i.e. start > end)
	rvs = coords$Orientation == 'Rev'
	coords$cis.start[rvs] = coords$Start[rvs] + upstream
	coords$cis.end[rvs] = coords$Start[rvs] - downstream
	regions=coords[,names(coords) %in% c('canonical_Name','where','cis.start','cis.end')]
	names(regions)=c('name','seq','start','end')
}

annotate.promoters =
	function(posns,regns)
{
	return( sapply(posns,
		function(x){
			regns.contain = regns$start <= x & regns$end >= x
			if(length(which(regns.contain))==0){regns.contain=NA}
			else{regns.contain=paste(regns$name[regns.contain],sep=',',collapse=',')}
			return(regns.contain)
		}
	))
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
	h$promoters = annotate.promoters((h[,4]+h[,5])/2,regions)

	# write file
	write.table(h,'hitsflat.tsv',quote=F,row.names=F,sep='\t')

} else {
	path = 'hitsflat.tsv'
	# this file should reflect everything done above
	h = read.delim(path)
}

# ANALYSIS/PLOTTING/OUTPUT
quit('no')
png=1

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

if(png){png('motif.hits.in.promoters.png',w=640,h=640)}
par(mar=c(10,10,10,6))
barplot(counts2,las=2,legend=c('not promoter','promoter region'),col=c('darkblue','red'),ylab='number of motif hits',main='rosetta motif hits in Halo promoters')
if(png){dev.off()}

p.enrich=counts2[2,]/(counts[1,]+counts2[2,])
p.enrich=sort(p.enrich,decreas=T)

if(png){png('motif.hits.in.promoters.fraction.png',w=640,h=640)}
par(mar=c(10,10,10,6))
barplot(p.enrich,las=2,ylab='fraction of motif hits',main='fraction of rosetta motif hits in Halo promoter regions\n(fraction of all genome hits per TF)')
if(png){dev.off()}

if(png){png('motif.hits.in.promoters.p.values.png',w=640,h=640)}
par(mar=c(10,10,10,6))
cols=c('darkblue','red')
boxplot(log(h$hit_p.value,10) ~ h$in.region * h$motif,las=2,col=cols,outcol=cols,outpch=21,legend=c('not promoter','promoter'), horizontal=T,main='p-values for motifs in promoters vs. not in promoters')
if(png){dev.off()}
