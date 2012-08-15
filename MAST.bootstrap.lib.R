#!/usr/bin/env Rscript

# MAST target sequence (e.g. genome) with a starting MEME motif,
#1 find matches (to do: optionally filtered to certain regions, e.g. promoters or ChIP overlaps?),
#2 build new motif based on actual matched sequences
#3 re-MAST with new [mixture] motif (bootstrap)
# repeat 1-3 until satisfied

load.sequence =
	# convenience function for loading fasta sequence into simple list format
	function(file)
{
	library(Biostrings)
	seqs=readFASTA(file,strip.descs=T)
	# convert from awkward list of pairs structure to more intuitive list structure
	seqlist=list()
	seqlist=sapply(seqs,function(x){x$seq})
	names(seqlist)=sapply(seqs,function(x){x$desc})
	return(seqlist)
}

get.hitseqs =
	function(hits,seqs)
{
	sapply(1:nrow(hits),
		function(i){
			contig = hits[i,1]
			strand = hits[i,2]
			start = hits[i,3]
			end = hits[i,4]
			seq = substring(seqs[[contig]],start,end)
#			if(strand==-1){ revcompDNA(seq) }
			if(strand==-1){ seq = as.character(reverseComplement(DNAString(seq))) }
			return(seq)
		}
	)
}

seqs.to.count.matrix =
	function(seqs)
{
	ltrs = sort(unique(unlist(strsplit(seqs,''))))
	maxl = max(sapply(seqs,nchar))
	cmat = matrix(0,nrow=maxl,ncol=length(ltrs))
	colnames(cmat) = ltrs
	# [sapply didn't want to increment in cmat here...]
	for(seq in seqs){
		for(i in 1:nchar(seq)){
			ltr = substring(seq,i,i)
			cmat[i,which(ltrs==ltr)] = cmat[i,which(ltrs==ltr)] + 1
		}
	}
	return(cmat)
}

# open starting MEME file, load prob matrix and display
# assumes implicit MEME formatting and 4-column prob matrix
MEME.pmat.from.file =
	function(fname)
{
	lines = readLines(fname)
	begin = FALSE
	pmat = c()
	for(line in lines){
		words = unlist(strsplit(line,'[ \t]'))
#		print(words)
		if( length(words)<1 | all(words==' ') | all(words=='\t') ){
#			cat('empty line\n')
			if(begin){return(matrix(pmat,ncol=4,byrow=T))}
			else next
		}
		if( words[1] == 'letter-probability' ){
#			cat('start\n')
			begin = TRUE
			next
		}
		if(begin){
			pmat = c(pmat,as.numeric(words))
		}
	}
}

write.MEME.motif =
	function(prob.matrix,bg.probs,mname,nhits)
{
	# MEME PWM
	weight.matrix = t(apply(prob.matrix,1,function(x){log(x/as.numeric(bg.probs),2)}))

	lines = c( "MEME version 3.0", "", "ALPHABET= ACGT", "", "strands: + -", "","Background letter frequencies (from dataset with add-one prior applied):" )
	lines = c( lines, paste( names(unlist(bg.probs)), sprintf( "%.3f", unlist( bg.probs )), collapse=" " ) )
	lines = c( lines, "", sprintf( "MOTIF %s", mname ), sprintf( "log-odds matrix: alength= %d w= %d", length(bg.probs), nrow(prob.matrix)))
	lines = c( lines, apply( weight.matrix, 1, function(i) sprintf( "%5.3f %5.3f %5.3f %5.3f", i[1],i[2],i[3],i[4] )))
	lines = c( lines, "", sprintf( "letter-probability matrix: alength= %d w= %d nsites= %d E= %.3e", length(bg.probs), nrow(prob.matrix), nhits, 0 ))
	lines = c( lines, apply( prob.matrix, 1, function(i) sprintf( "%5.3f %5.3f %5.3f %5.3f", i[1],i[2],i[3],i[4] )))

	fname = paste(mname,'MEME',sep='.')
	cat( lines, sep="\n", file=fname )
	return(fname)
}

positions.in.regions =
	function(positions,p1p2)
{
	return(sapply(p1p2,
		function(x){
			x >= p1p2[,1] & x <= p1p2[,2]
		}
	))
}

default_bgprobs =
	function()
{
	bgprobs = list(
		'A'=0.25,
		'C'=0.25,
		'G'=0.25,
		'T'=0.25
	)
	return(bgprobs)
}

load_bgprobs =
	function(file)
{
	bgprobs = list()
	for( line in readLines(file) ) {
		fields = unlist(strsplit(line,'[ \t]'))
		bgprobs[fields[1]] = as.numeric(fields[2])
	}
	return(bgprobs)
}

load_regions_from_genecoords =
	function(file,window=c(300,100))
{
# load regions in order to filter genome hits
# regions is a 2-column matrix of start and end positions, with optional rownames
	upstream=window[1]
	downstream=window[2]
	coords=read.delim(file,as.is=T)
	coords$Start = as.numeric(coords$Start)
	coords$cis.start = 0
	coords$cis.end = 0
	fwd = coords$Orientation == 'For'
	coords$cis.start[fwd] = coords$Start[fwd] - upstream
	coords$cis.end[fwd] = coords$Start[fwd] + downstream
	# in halo gene coords, 'Start' (generally) refers to coordinate start (i.e. start > end: not 5' for reverse-strand genes)
	rvs = coords$Orientation == 'Rev'
	# is this right? might need to be looking at coords$End
	coords$cis.start[rvs] = coords$Start[rvs] + upstream
	coords$cis.end[rvs] = coords$Start[rvs] - downstream
	regions=coords[, c('canonical_Name','where','cis.start','cis.end')]
	names(regions)=c('name','seq','start','end')
#	cat(nrow(regions), ' regions\n')
	print( head(regions) )
	return(regions)
}

load_regions_from_TSS =
	function(file,window=c(300,100))
{
# load regions in order to filter genome hits
# regions to return is a 2-column matrix of start and end positions, with optional rownames
	upstream=window[1]
	downstream=window[2]
	TSS=read.delim(file,as.is=T)
	TSS = TSS[ !is.na(TSS$TSS), ]
	TSS$TSS = as.integer( as.numeric(TSS$TSS) )
	TSS$window.start = 0
	TSS$window.end = 0
	fwd = TSS$Strand == 'For'
	TSS$window.start[fwd] = TSS$TSS[fwd] - upstream
	TSS$window.end[fwd] = TSS$TSS[fwd] + downstream
	# in halo gene TSS, 'Start' refers transcriptional start, even for reverse-stranded genes
	rvs = TSS$Strand == 'Rev'
	TSS$window.start[rvs] = TSS$TSS[rvs] + upstream
	TSS$window.end[rvs] = TSS$TSS[rvs] - downstream
	regions=TSS[, c('Name','Chr','window.start','window.end')]
	names(regions)=c('name','seq','start','end')
#	cat(nrow(regions), ' regions\n')
	print( head(regions) )
	return(regions)
}

regions_for_genes_file =
	function(file,regions)
{
	genes = readLines(file)
	return( regions_for_genes(genes,regions) )
}

regions_for_genes =
	function(genes,regions)
{
#	cat('filtering ', nrow(regions), ' regions down to ', length(genes), ' genes\n')
	regions = regions[ regions$name %in% genes, ]
#	cat('there are now ', nrow(regions), ' regions\n')
	return(regions)
}

filter.hits.by.regions =
	function(hits,regions)
{
	sapply(1:nrow(hits),
		function(h){
			seq.match = hits[h,1] == regions$seq
			pos = (hits[h,3] + hits[h,4]) / 2
			after.start = pos >= regions$start
			before.end = pos <= regions$end
			return( TRUE %in% (seq.match & after.start & before.end) )
		}
	)
}

motif_bootstrap =
	function(motfile,seqfile,bgfile,regionsfile,genelist=NULL,niter=10,maxhits=500,mixture=0.25,pseudocounts=1,mt=0.0001,images=T,prefix='BS.')
{
	logf = paste(prefix,'log',sep='/')
	cat(motfile,seqfile,bgfile,regionsfile,niter,maxhits,mixture,pseudocounts,mt,images,'\n',file=logf)

	mast.args = sprintf('-hit_list -mt %g',mt)
	bg.probs = default_bgprobs()
	if(!is.null(bgfile)){
		mast.args = paste(mast.args,'-bfile',bgfile)
		bg.probs = load_bgprobs(bgfile)
	}

	# NOTE: 0001 will be the starting search motif, 0002 will be the first bootstrap motif
	if(images){
		pdf(paste(prefix,'motif.%04d.pdf',sep=''),width=10,height=4,onefile=F)
	}

	start.pmat = MEME.pmat.from.file(motfile)
	library(seqLogo)
	seqLogo(t(start.pmat))

	if(!images) dev.new()

	last.prob.matrix = NULL

	#regions = load_regions_from_genecoords(regionsfile)
	regions = load_regions_from_TSS(regionsfile)

	if (!is.null(genelist)){
		cat('filtering ', nrow(regions), ' regions down to ', length(genelist), ' genes\n', file=logf, append=T)
		regions = regions_for_genes(genelist,regions)
	}
	cat('there are ', nrow(regions), ' regions\n', file=logf, append=T)

	# TO DO: compute true background in regions being searched?

	for( iter in 1:niter ) {
		cat('iter ',iter,'\n',sep='')
		cat('MASTing motif from file "',motfile,'" in sequence from fasta file "',seqfile,'"\n',sep='')
		hitfile = sprintf('%shitfile.%04d',prefix,iter,sep='')
		mast.cmd = paste('mast',motfile,seqfile,mast.args,'>',hitfile)
		cat(mast.cmd,'\n')
		system(mast.cmd)
		hits = read.table(hitfile,as.is=T)

		# filter by regions (if present)
		if(!is.null(regions)) hits = hits[ filter.hits.by.regions(hits,regions), ]
		cat(nrow(hits),' hits in regions\n',file=logf,append=T)

		if (nrow(hits) == 0) {
			cat('no hits in regions--failed to find/adapt motif!\n',file=logf,append=T)
			break
		}

		# sort by pval and cull
		hits = hits[ order(hits[,6])[ 1:min(nrow(hits),maxhits) ] , ]

		# load searched sequence(s)
		seqs = load.sequence(seqfile)

		# get sequences corresponding to matches
		hitseqs = get.hitseqs(hits,seqs)
		cat(paste(head(hitseqs),'\n',sep=''),'...\n',sep='')

		count.matrix = seqs.to.count.matrix(hitseqs)
		# add pseudocounts for mathematical stability
		count.matrix[ count.matrix == 0 ] = pseudocounts
		print(count.matrix)

		# probability matrix
		prob.matrix = count.matrix / apply(count.matrix,1,sum)
		if ( !all(prob.matrix==0.25) ) seqLogo(t(prob.matrix))
		if ( is.null(last.prob.matrix) ) {
		last.prob.matrix = prob.matrix
		} else if ( all(prob.matrix == last.prob.matrix) ) {
			cat('bootstrap PPM converged!\n',file=logf,append=T)
			break
		} else {
			# mix the match-generated motif with the search motif
			# how? lets just try weighted average probabilities...
			prob.matrix = (1-mixture) * prob.matrix + mixture * last.prob.matrix
		}

		mname = sprintf('%s%s.%04d',prefix,seqfile,iter)
		# write new bootstrap motif to file, updating 'motfile' in the process so that it is used next
		motfile = write.MEME.motif(prob.matrix,bg.probs,mname,length(hitseqs))
	}

	if (images) dev.off()
}
