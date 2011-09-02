annote.promoters =
	function(positions,upstream=400,downstream=50,gene.coords='halo.gene.coords.tsv')
{
	cat('annotating positions in cis regions, window [',-1*upstream,',',downstream,']\n')
	# halo gene coords
	gc=read.delim(gene.coords)

	# better: halo TSS where known
	# tss = read.delim('')

	gc$Start = as.numeric(gc$Start)
	gc$cis.start = 0
	gc$cis.end = 0
	fwd = gc$Orientation == 'For'
	gc$cis.start[fwd] = gc$Start[fwd] - upstream
	gc$cis.end[fwd] = gc$Start[fwd] + downstream
	# in halo gene coords, 'Start' (generally) refers to actual start (i.e. start > end)
	rvs = gc$Orientation == 'Rev'
	gc$cis.start[rvs] = gc$Start[rvs] + upstream
	gc$cis.end[rvs] = gc$Start[rvs] - downstream

	positions$in.region=NA

	sequences = levels(factor(gc$where))
	for(s in sequences){
		cat('sequence',s,'\n')
		seq.genes = which(gc$where==s)
		positions$in.region[positions$sequence_name==s] =
			sapply(positions$hit_start[positions$sequence_name==s], function(pos){ TRUE %in%
				sapply(seq.genes, function(g){pos > gc$cis.start[g] & pos < gc$cis.end[g]})
			})
	}
	return(positions)
}
