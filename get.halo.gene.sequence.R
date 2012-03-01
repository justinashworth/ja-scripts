# for load.sequence
source('~/scripts/isb-scripts/misc.R')

get.halo.gene.seq =
	function(
		vngids,
		seqrange=c(0,NA),
		coords='~/baliga/halo/genome/halo.gene.coords.tsv',
		gseq='~/baliga/halo/genome/gene_sequence.fa'
	)
{
	if(!exists('halo.gene.coords')){
		halo.gene.coords=read.delim(coords)
		# assign to global environment so we don't have to keep loading the file
		assign('halo.gene.coords',halo.gene.coords,envir=.GlobalEnv)
	}
	if(!exists('halo.gseq')){
		require(Biostrings)
		halo.gseq=load.sequence(gseq)
		assign('halo.gseq',halo.gseq,envir=.GlobalEnv)
	}
	indices = sapply( vngids, function(x) { grep(x,halo.gene.coords$canonical_Name) } )
	contigs = halo.gene.coords$where[indices]
	# note about stops and starts in halo.gene.coords file: Start > Stop for Rev strand genes
	starts = halo.gene.coords$Start[indices]
	stops  = halo.gene.coords$Stop [indices]
	fwds = halo.gene.coords$Orientation[indices] == 'For'
	rvss = halo.gene.coords$Orientation[indices] == 'Rev'
	# adjust for custom sequence range (relative to start)
	starts[fwds] = starts[fwds] + seqrange[1]
	starts[rvss] = starts[rvss] - seqrange[1]
	if ( ! is.na(seqrange[2]) ) {
		stops = halo.gene.coords$Start[indices]
		stops[fwds] = stops[fwds] + seqrange[2]
		stops[rvss] = stops[rvss] - seqrange[2]
	}
	# flip Rev starts and stops for substring
	stops2 = stops
	stops[rvss] = starts[rvss]
	starts[rvss] = stops2[rvss]

	seqs = sapply( 1:length(indices), function(i) { substring(halo.gseq[[ contigs[i] ]], starts[i] ,stops[i] ) } )
	if ( any(rvss) ) seqs[rvss] = sapply( seqs[rvss], function(x) { as.character( reverseComplement( DNAString(x) ) ) } )
	names(seqs) = vngids
	return(seqs)
}
