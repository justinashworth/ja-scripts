# search target sequence (e.g. genome) with a starting MEME motif,
#1 find matches
#2 build new motif based on actual matched sequences
#3 re-search with new [mixture] motif
# repeat 1-3 until converged

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

rvscomp =
	function(x)
{
	return( as.character(reverseComplement(DNAString(x))) )
}

get.hitseqs =
	function(hits,seqs,flank=c(0,0),revcomp=T)
{
	sapply(1:nrow(hits),
		function(hit){
			contig = hits[hit,'seq']
			strand = hits[hit,'strand']
			start = hits[hit,'start']
			end = hits[hit,'end']
			if(!contig %in% names(seqs)){
				cat('in get.hitseqs:',contig,' not in sequences!\n')
			}
			seq = substring(seqs[[contig]], start+flank[1], end+flank[2])
			# sometimes, decision to transform is external (revcomp==F)
			if(strand==-1 & revcomp){seq = rvscomp(seq)}
			return(seq)
		}
	)
}

seqs.to.count.matrix =
	function(seqs,dna=T)
{
	ltrs = c('A','C','G','T')
	if(!dna) ltrs = sort(unique(unlist(strsplit(seqs,''))))

	maxlen = max(sapply(seqs,nchar))
	cmat = matrix(0,nrow=maxlen,ncol=length(ltrs))
	colnames(cmat) = ltrs
	# [sapply didn't want to increment in cmat here...]
	for(seq in seqs){
		for(char in 1:nchar(seq)){
			ltr = substring(seq,char,char)
			if(!ltr %in% ltrs) next
			cmat[char,ltr] = cmat[char,ltr] + 1
		}
	}
	return(cmat)
}

weighted.count.matrix =
	function(hits,weights,sequence)
{
	hitseqs = get.hitseqs(hits,sequence)
	names(hitseqs) = hits$regionnames
	cat(paste(head(hitseqs),'\n',sep=''),'...\n',sep='')

	ltrs = c('A','C','G','T')
	maxlen = max(sapply(hitseqs,nchar))
	cmat = matrix(0,nrow=maxlen,ncol=length(ltrs))
	colnames(cmat) = ltrs

	for(gene in names(hitseqs)){
		seq = hitseqs[[gene]]
		weight = 0
		if(gene %in% names(weights)){
			weight = weights[[gene]]
		}
		for(char in 1:nchar(seq)){
			ltr = substring(seq,char,char)
			if(!ltr %in% ltrs) next
			cmat[char,ltr] = cmat[char,ltr] + weight
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
	function(prob.matrix,bg.probs,mname,nhits=1)
{
	# MEME PWM
	weight.matrix = t(apply(prob.matrix,1,function(x){log(x/as.numeric(bg.probs),2)}))

	lines = c( "MEME version 3.0", "", "ALPHABET= ACGT", "", "strands: + -", "","Background letter frequencies (from dataset with add-one prior applied):" )
	lines = c( lines, paste( names(unlist(bg.probs)), sprintf( "%.4f", unlist( bg.probs )), collapse=" " ) )
	lines = c( lines, "", sprintf( "MOTIF %s", mname ), sprintf( "log-odds matrix: alength= %d w= %d", length(bg.probs), nrow(prob.matrix)))
	lines = c( lines, apply( weight.matrix, 1, function(i) sprintf( "%5.3f %5.3f %5.3f %5.3f", i[1],i[2],i[3],i[4] )))
	lines = c( lines, "", sprintf( "letter-probability matrix: alength= %d w= %d nsites= %d E= %.3e", length(bg.probs), nrow(prob.matrix), nhits, 0 ))
	lines = c( lines, apply( prob.matrix, 1, function(i) sprintf( "%.4f %.4f %.4f %.4f", i[1],i[2],i[3],i[4] )))

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
	function(file,window=c(250,50))
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
	coords$cis.start[fwd] = as.integer(coords$Start[fwd] - upstream)
	coords$cis.end[fwd] = as.integer(coords$Start[fwd] + downstream)
	# in halo gene coords, 'Start' is the start of transcription, and thus for reverse-strand genes is > 'Stop'
	# for our purposes however, we will define region start < end, even for rvs genes (this makes region checks simpler)
	rvs = coords$Orientation == 'Rev'
	coords$cis.start[rvs] = as.integer(coords$Start[rvs] - downstream)
	coords$cis.end[rvs] = as.integer(coords$Start[rvs] + upstream)
	regions=coords[, c('canonical_Name','where','cis.start','cis.end','Orientation','Start','Stop')]
	names(regions)=c('name','seq','start','end','dir','tss','tts')
	#regions=coords[, c('canonical_Name','where','cis.start','cis.end','Orientation')]
	#names(regions)=c('name','seq','start','end','dir')
#	cat(nrow(regions), ' regions\n')
#	print( head(regions) )
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
	regions=TSS[, c('Name','Chr','window.start','window.end','Strand','TSS','TTS')]
	names(regions)=c('name','seq','start','end','dir','tss','tts')
#	cat(nrow(regions), ' regions\n')
	print( head(regions) )
	return(regions)
}

seqs_for_regions =
	function(regs,seqs='halo.genome.fa')
{
	require(Biostrings)
	gs = readDNAStringSet(seqs)
	apply(regs, 1, function(x){
		cat(x[['name']],'\n')
		start = as.integer(x[['start']])
		end = as.integer(x[['end']])
		cat(start,end,'\n')
		if(start<1) start = 1
		ll=length(gs[[x[['seq']]]])
		if(end>ll) end = ll
		cat(start,end,'\n')
		ss=gs[[ x[['seq']] ]][ start:end ]
		if(x[['dir']]=='Rev'){ ss=reverseComplement(ss) }
		return(ss)
	})
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

hits.in.regions =
	function(hits,regions)
{
	# filter/annotate hits in regions,
	# expand/flatten for hits matching multiple regions as necessary
	hitindices = c()
	regionnames = c()

	if(!all(regions$end >= regions$start)){
		cat('PROBLEM: region starts > region ends\n')
		log(0)
	}

	for(hit in 1:nrow(hits)){
		seq.match = as.character(hits[hit,'seq']) == as.character(regions$seq)
		pos = (hits[hit,'start'] + hits[hit,'end']) / 2
		after.start = pos >= regions$start
		before.end = pos <= regions$end
		matched.regions = na.exclude( regions$name[ which(seq.match & after.start & before.end) ] )
		if(length(matched.regions) < 1) next
		for(match in 1:length(matched.regions)){
			hitindices = c(hitindices, hit)
			regionnames = c(regionnames, as.character(matched.regions[match]))
		}
	}
	hitsinregions = unique( cbind( hits[hitindices,], regionnames) )
	return(hitsinregions)
}

motif_refine = function(
	motfile,
	seqfile,
	bgfile,
	regionsfile,
	genelist=NULL,
	niter=100,
	maxhits=500,
	mixture_start=0.75,
	mixture_end=0.25,
	finish_iters=10,
	pseudocounts=1,
	mt_start=0.0005,
	mt_end=0.0001,
	images=T,
	prefix='MR.',
	geneweights=NULL
){
	result = list(motif=NULL, hits=NULL, target.genes=NULL, target.genes.matched=NULL, all.genes=NULL, all.genes.matched=NULL, enrichment.pval=NULL)
	logf = paste(prefix,'log',sep='/')
	cat(motfile,seqfile,bgfile,regionsfile,niter,maxhits,mixture_start,mixture_end,finish_iters,pseudocounts,mt_start,mt_end,images,'\n',file=logf)

	# load searched sequence(s)
	search_sequence = load.sequence(seqfile)
	cat('loaded sequence file with contigs:', names(search_sequence), '\n')

	#mast.args = sprintf('-hit_list -mt %g',mt)
	#mast.args = '-hit_list'
	fimo.args = '--text'
	bg.probs = default_bgprobs()
	if(!is.null(bgfile)){
		#mast.args = paste(mast.args,'-bfile',bgfile)
		fimo.args = paste(fimo.args,'-bgfile',bgfile)
		bg.probs = load_bgprobs(bgfile)
	}

	start.pmat = MEME.pmat.from.file(motfile)
	last.prob.matrix = start.pmat
	# rewrite starting motif to ensure faithfulness of following function
	mname = sprintf('%s%s.%04d',prefix,seqfile,0)
	write.MEME.motif(start.pmat,bg.probs,mname)

	library(seqLogo)

	if(images){ pdf( sprintf('%smotif.%04d.pdf',prefix,0), width=10,height=4) }

	seqLogo(t(start.pmat))

	if(images) dev.off()
	else dev.new()

	allregions = load_regions_from_genecoords(regionsfile)
	n.allregions = length( unique(allregions$name) )
	#allregions = load_regions_from_TSS(regionsfile)

	regions = allregions

	# filter regions down to genes of interest
	if (!is.null(genelist)){
		genelist = unique(genelist)
		cat('filtering ', nrow(regions), ' regions down to ', length(genelist), ' genes\n', file=logf, append=T)
		regions = regions_for_genes(genelist,regions)
	}
	n.regions = length( unique(regions$name) )
	cat('there are ', n.regions, ' regions\n', file=logf, append=T)

	for( iter in 1:(niter+finish_iters) ) {
		if(iter <= niter){
			mixture = mixture_start - (mixture_start - mixture_end) * (iter-1)/(niter-1)
			mt = mt_start - (mt_start - mt_end) * (iter-1)/(niter-1)
		} else {
			# finish_iters (allow motifs to settle with fixed final mixture factor and motif threshold)
			mixture = mixture_end
			mt = mt_end
		}
		cat('iter',iter,'mixture',mixture,'mt',mt,'\n')
		#cat('MASTing motif from file "',motfile,'" in sequence from fasta file "',seqfile,'"\n',sep='')
		cat('FIMO with motif from file "',motfile,'" in sequence from fasta file "',seqfile,'"\n',sep='')
		#masthits = sprintf('%smasthits.%04d',prefix,iter,sep='')
		motifhits = sprintf('%smotifhits.%04d',prefix,iter,sep='')
		#mast.cmd = sprintf('mast_4.9.0 %s %s %s -mt %g > %s',motfile,seqfile,mast.args,mt,masthits)
		search.cmd = sprintf('fimo_4.9.0 %s --thresh %g %s %s > %s',fimo.args,mt,motfile,seqfile,motifhits)
		cat(search.cmd,'\n')
		system(search.cmd)
		#hits = read.table(motifhits,as.is=T)
		#names(hits) = c('seq','strand','start','end','score','pval')
		hits = read.delim(motifhits,as.is=T,comment.char='#',header=F)
		names(hits)=c('motif', 'seq', 'start', 'end', 'strand', 'score', 'pval', 'q-value', 'matched_sequence')
		# filter hits by p-value
		hits = hits[ hits$pval <= mt, ]
		#print(head(hits,10))

		# annotate all hits (genomewide)
		# warning: slow. Avoid if not necessary
		allhits = hits.in.regions(hits,allregions)
		allhitsfile = sprintf('%sallhits.%04d',prefix,iter,sep='')
		#tsv(allhits,allhitsfile)

		allhitgenesfile = sprintf('%sallgeneshit.%04d',prefix,iter,sep='')
		#cat( unique(as.character(allhits$regionnames)), file=allhitgenesfile)

		# TODO get number of all gene regions containing hits
		n.allregions.hit = length( unique(allhits$regionnames) )

		# filter hits by regions (if present)
		if(!is.null(regions)) {
			cat('fetching hits in regions\n')
			hits = hits.in.regions(hits,regions)

			tgtgeneshitfile = sprintf('%stargetgeneshit.%04d',prefix,iter,sep='')
			#cat( unique(as.character(hits$regionnames)), file=tgtgeneshitfile)

			n.regions.hit = length( unique(hits$regionnames) )
			#regionhits = sprintf('%sregionhits.%04d',prefix,iter,sep='')
			#tsv(hits,regionhits)

			# compute a hypergeometric p for hits in gene vs. hits in genome
			cat(n.regions.hit, 'regions hit in', n.regions, 'target regions with', n.allregions.hit, 'of all', n.allregions, 'total regions hit\n')
			p.regionhits = phyper(n.regions.hit, n.allregions.hit, n.allregions-n.allregions.hit, n.regions, lower.tail=F)
			result$enrichment.pval = p.regionhits
			cat('hypergeometric p-val:', p.regionhits, '\n')

			cat(n.regions.hit,'/',n.regions,' hits in target gene upstream regions with ', n.allregions.hit,'/',n.allregions, ' of all regions hit', ' pval=', p.regionhits, '\n',sep='',file=logf,append=T)
		}
		result$hits = hits
		result$target.genes = genelist
		result$target.genes.matched = unique(as.character(hits$regionnames))
		result$all.genes = unique(allregions$name)
		result$all.genes.matched = unique(as.character(allhits$regionnames))

		regions = regions[ !is.na(regions$seq), ]

		if (nrow(hits) == 0) {
			cat('no hits in regions--failed to find/adapt motif!\n',file=logf,append=T)
			break
		}

		# sort by pval and cull
		hits = hits[ order(hits[,'pval'])[ 1:min(nrow(hits),maxhits) ] , ]

		# write sorted, culled hits
		besthits = sprintf('%sbesthits.%04d',prefix,iter, sep='')
		#tsv(hits,besthits)
		print(head(hits,10))

		# get sequences corresponding to matches
		hitseqs = get.hitseqs(hits,search_sequence)
		cat(paste(head(hitseqs),'\n',sep=''),'...\n',sep='')

		if(!is.null(geneweights)){
			# make a pseudo-count matrix weighted by orthogonal gene scores (e.g. ChIP strength/certainty)
			count.matrix = weighted.count.matrix(hits,geneweights,search_sequence)
			cat('Weighted count matrix (raw):\n')
			print(count.matrix)
		} else {
			count.matrix = seqs.to.count.matrix(hitseqs)
			cat('Count matrix (raw):\n')
			print(count.matrix)
		}

		## add pseudocounts for mathematical stability
		count.matrix = count.matrix + pseudocounts
		#cat('Count matrix (with pseudocounts):\n')
		#print(count.matrix)

		# probability matrix
		prob.matrix = count.matrix / apply(count.matrix,1,sum)
		if ( all(prob.matrix == last.prob.matrix) ) {
			cat('PPM converged!\n',file=logf,append=T)
			break
		} else {
			# mix the match-generated motif with the search motif
			# weighted average probabilities
			prob.matrix = (1-mixture) * prob.matrix + mixture * last.prob.matrix
		}

		if ( !all(prob.matrix==0.25) ) {
			if(images){ pdf(sprintf('%smotif.%04d.pdf',prefix,iter),width=10,height=4) }
			seqLogo(t(prob.matrix))
			if(images) dev.off()
			else dev.new()
		}

		mname = sprintf('%s%s.%04d',prefix,seqfile,iter)
		# write new motif to file, updating 'motfile' in the process so that it is used next
		motfile = write.MEME.motif(prob.matrix,bg.probs,mname,length(hitseqs))
		result$motif = prob.matrix

		# plot best aligned hit sites in sequences of interest
		if(images & (iter==(niter+finish_iters) | iter==1)){

			hit.centric = F
			tss.centric = T

			if(hit.centric){

				pdf( sprintf('%salignedhits.%04d.pdf', prefix, iter), onefile=F, useDingbats=F, font='mono', pointsize=8, height=8, width=28)
				par(family='mono')
				nhits = min(50, nrow(hits))
				xbound = c(-100,300)
				plot(0, 0, type='n', xlim=c(xbound[1]-10,xbound[2]), ylim=c(-1,nhits+1), ylab='', xlab='', yaxt='n', main=mname, xaxs='i', yaxs='i')

				for(hit in 1:nhits){

					genelab = hits[hit,'regionnames']
					region = which(regions$name==genelab)[1]
					#print(regions[region,])
					tss = as.integer( regions[region,'tss'] )
					tts = as.integer( regions[region,'tts'] )
					rvs = regions[region,'dir'] == 'Rvs'
					if(is.na(tts)){
						#...it happens
						if(rvs)tts = tss-500
						else tts = tss+500
					}
					#cat(genelab,tss,tts,rvs,'\n')

					# hit sequence
					hitseq = get.hitseqs(hits[hit,],search_sequence,revcomp=F)
					if(rvs){hitseq = rvscomp(hitseq)}
#					print(hits[hit,])
#					print(hitseq)

					# region sequence
					flank = c(-200,200)
					bgseq = get.hitseqs(hits[hit,],search_sequence,flank=flank,revcomp=F)
					if(rvs){bgseq = rvscomp(bgseq)}

					for(char in 1:nchar(bgseq)){
						xpos = char-1+flank[1]
						if(xpos < xbound[1] | xpos > xbound[2]) next
						if(char>=1+flank[1] & char<=nchar(hitseq)+flank[1])next
						letter = substring(bgseq,char,char)
						text(xpos, nhits-hit, letter, col='black',font=1)
					}

					# gene label and TSS
					text(-60, nhits-hit, genelab)
					hitstart = hits[hit,'start']
					if(!rvs){
						offset = tss-hitstart
						tsslen = tts-tss
						maxchar = min(tsslen,200)
						tssseq = get.hitseqs(hits[hit,],search_sequence,flank=c(offset,offset+maxchar),revcomp=F)
						for(char in 1:nchar(tssseq)){
							xpos = char-1+offset
							if(xpos < xbound[1] | xpos > xbound[2]) next

							#if(char>=1-flank[1] & char<=nchar(hitseq)-flank[1])next
							letter = substring(tssseq,char,char)
							text(xpos, nhits-hit, letter, col='orange',font=2)
						}
					}
					# hit seq (top layer)
					for(char in 1:nchar(hitseq)){
						letter = substring(hitseq,char,char)
						text(char-1, nhits-hit, letter, col='blue',font=2)
					}
				}
				dev.off()
			}

			if(tss.centric){
				pdf( sprintf('%salignedtss.%04d.pdf', prefix, iter), onefile=F, useDingbats=F, font='mono', pointsize=8, height=8, width=28)
				par(family='mono')

				nhits = min(50, nrow(hits))
				topgenes = regions[ regions$name %in% hits$regionnames[1:nhits], ]
				ngenes.plot = nrow(topgenes)

				xbound = c(-300,100)
				plot(0, 0, type='n', xlim=c(xbound[1]-10,xbound[2]), ylim=c(-5,ngenes.plot+5), ylab='', xlab='', yaxt='n', main=mname)

				# plot array of these sequences, aligned by the TSS
				for(gene in 1:ngenes.plot){
					name = topgenes[gene,'name']
					text(xbound[1]-10, ngenes.plot-gene, name)

					tss = as.integer( topgenes[gene,'tss'] )
					tts = as.integer( topgenes[gene,'tts'] )
					rvs = topgenes[gene,'dir'] == 'Rvs'
					if(is.na(tts)){
						#...it happens
						if(rvs)tts = tss-500
						else tts = tss+500
					}

					#print(topgenes[gene,])
					contig = topgenes[gene,'seq']
					upstream = character()
					downstream = character()

					if(!contig %in% names(search_sequence)){
						cat('in tss plotting:',contig,' not in sequences!\n')
					}
					if(rvs){
						upstream = rvscomp( substring(search_sequence[[contig]], tss+1, tss-xbound[1]))
						downstream = rvscomp( substring(search_sequence[[contig]], tss-xbound[2], tss))
					}else{
						upstream = substring(search_sequence[[contig]], tss+xbound[1], tss-1)
						downstream = substring(search_sequence[[contig]], tss, tss+xbound[2])
					}
					for(char in 1:nchar(upstream)){
						xpos = char+xbound[1]-1
						if(xpos < xbound[1] | xpos > xbound[2]) next
						letter = substring(upstream, char, char)
						text(xpos, ngenes.plot-gene, letter, col='black')
					}

					for(char in 1:nchar(downstream)){
						xpos = char-1
						if(xpos < xbound[1] | xpos > xbound[2]) next
						letter = substring(downstream, char, char)
						text(xpos, ngenes.plot-gene, letter, col='orange', font=2)
					}

					# superimpose hit sites for this gene
					genehits = hits[ hits$regionnames==name, ]
					for(hit in 1:nrow(genehits)){
						hitseq = get.hitseqs(genehits[hit,],search_sequence,revcomp=F)
						if(rvs) hitseq = rvscomp(hitseq)

						for(char in 1:nchar(hitseq)){
							xpos = char+genehits[hit,'start']-tss-1
							if(rvs) xpos = char+genehits[hit,'end']-tss-1
							letter = substring(hitseq, char, char)
							text(xpos, ngenes.plot-gene, letter, col='blue', font=2)
						}

					}
				}
				dev.off()
			}

		}

	}
	return(result)

}
