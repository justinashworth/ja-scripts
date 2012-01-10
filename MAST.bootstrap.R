#!/usr/bin/env Rscript

# MAST target sequence (e.g. genome) with a starting MEME motif,
#1 find matches (to do: optionally filtered to certain regions, e.g. promoters or ChIP overlaps?),
#2 build new motif based on actual matched sequences
#3 re-MAST with new [mixture] motif (bootstrap)
# repeat 1-3 until satisfied

# example
#./MAST.bootstrap.R -m T.vol.FL3.MEME -s halo.genome.fa -n 100 -x 500 -b halo.bg.file -p 10 -g -X 0.25

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
	function(prob.matrix,bg.probs,mname)
{
	# MEME PWM
	weight.matrix = t(apply(prob.matrix,1,function(x){log(x/as.numeric(bg.probs),2)}))

	lines = c( "MEME version 3.0", "", "ALPHABET= ACGT", "", "strands: + -", "","Background letter frequencies (from dataset with add-one prior applied):" )
	lines = c( lines, paste( names(unlist(bg.probs)), sprintf( "%.3f", unlist( bg.probs )), collapse=" " )  )
	lines = c( lines, "", sprintf( "MOTIF %s", mname ), sprintf( "log-odds matrix: alength= %d w= %d", length(bg.probs), nrow(prob.matrix)))
	lines = c( lines, apply( weight.matrix, 1, function(i) sprintf( "%5.3f %5.3f %5.3f %5.3f", i[1],i[2],i[3],i[4] )))
	lines = c( lines, "", sprintf( "letter-probability matrix: alength= %d w= %d nsites= %d E= %.3e", length(bg.probs), nrow(prob.matrix), length(hitseqs), 0  ))
	lines = c( lines, apply( prob.matrix, 1, function(i) sprintf( "%5.3f %5.3f %5.3f %5.3f", i[1],i[2],i[3],i[4] )))

	fname = paste(mname,'MEME',sep='.')
	cat( lines, sep="\n", file=fname )
	return(fname)
}

# INIT
library(seqLogo)
library('getopt')
opt = getopt(matrix(c(
	'motf'         , 'm', 1, 'character',
	'seqf'         , 's', 1, 'character',
	'niter'        , 'n', 1, 'integer',
	'maxhits'      , 'x', 1, 'integer',
	'mixture'      , 'X', 1, 'double',
	'bgfile'       , 'b', 1, 'character',
	'regions'      , 'r', 1, 'character',
	'pseudocounts' , 'p', 1, 'integer',
	'images'       , 'g', 0, 'logical',
	'help'         , 'h', 0, 'logical'
),ncol=4,byrow=T))

if( is.null(opt$niter) ) opt$niter = 10
if( is.null(opt$maxhits) ) opt$maxhits = 500
if( is.null(opt$mixture) ) opt$mixture = 0.25
if( is.null(opt$pseudocounts) ) opt$pseudocounts = 1
if( is.null(opt$images) ) opt$images = FALSE

if( is.null(opt$motf) | is.null(opt$seqf) | !is.null(opt$help) ) {
	cat('args: -m <motif_file> -s <fasta_seq> [-[n|-niter] <%d>] [-[x|-maxhits] <%d>] [-[b|-bgfile] <bgfile>] [-[p|-pseudocounts] <%d>] [-[g|-images]]\n')
	opt$motf = 'T.vol.FL3.MEME'
	opt$seqf = 'halo.genome.fa'
	opt$bgfile = 'halo.bg.file'
	opt$regions = 'halo.gene.coords.tsv'
#	q('no')
}

if( opt$mixture < 0 | opt$mixture > 1 ){
	cat('invalid mixture ratio\n')
	q('no')
}

bg.probs = list()

if( is.null(opt$bgfile) ) {
	bg.probs = list(
		'A'=0.25,
		'C'=0.25,
		'G'=0.25,
		'T'=0.25
	)
} else {
	for( line in readLines(opt$bgfile) ) {
		fields = unlist(strsplit(line,'[ \t]'))
		bg.probs[fields[1]] = as.numeric(fields[2])
	}
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

print(opt)

cat('background probabilities:\n')
print(paste(unlist(bg.probs)),sep='')
cat('NOTE: currently stripping letter info in bg probs and assuming implicit A,C,G,T ordering!\n')

# load regions in order to filter genome hits
# regions is a 2-column matrix of start and end positions, with optional rownames
regions = NULL
if( !is.null(opt$regions) ){
	upstream=500
	downstream=100
#	coords=read.delim('halo.gene.coords.tsv')
	coords=read.delim(opt$regions,as.is=T)
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

mast.args = '-hit_list'
if( !is.null(opt$bgfile) ) { mast.args = paste(mast.args,'-bfile',opt$bgfile) }

# NOTE: 0001 will be the starting search motif, 0002 will be the first bootstrap motif
if(opt$images) pdf('BS.motif.%04d.pdf',width=10,height=4,onefile=F)

start.pmat = MEME.pmat.from.file(opt$motf)
seqLogo(t(start.pmat))

if(is.null(opt$images) | !opt$images) dev.new()

mfile = opt$motf
last.prob.matrix = NULL
for( iter in 1:opt$niter ) {
	cat('iter ',iter,'\n',sep='')
	cat('MASTing motif from file "',mfile,'" in sequence from fasta file "',opt$seqf,'"\n',sep='')
	hitfile = sprintf('BS.hitfile.%04d',iter,sep='')
	mast.cmd = paste('mast',mfile,opt$seqf,mast.args,'>',hitfile)
	cat(mast.cmd,'\n')
	system(mast.cmd)
	hits = read.table(hitfile,as.is=T)

	# filter by regions (if present)
	if(!is.null(regions)) hits = hits[ filter.hits.by.regions(hits,regions), ]
	cat(nrow(hits),' hits in regions\n')

	# sort by pval and cull
	hits = hits[ order(hits[,6])[ 1:min(nrow(hits),opt$maxhits) ] , ]

	# load searched sequence(s)
	seqs = load.sequence(opt$seqf)

	# get sequences corresponding to matches
	hitseqs = get.hitseqs(hits,seqs)
	cat(length(hitseqs),' matches:\n',paste(head(hitseqs),'\n',sep=''),'...\n',sep='')

	count.matrix = seqs.to.count.matrix(hitseqs)
	print(count.matrix)
	# add pseudocounts for mathematical stability
	count.matrix[ count.matrix == 0 ] = opt$pseudocounts

	# probability matrix
	prob.matrix = count.matrix / apply(count.matrix,1,sum)
	seqLogo(t(prob.matrix))
	if ( is.null(last.prob.matrix) ) {
	last.prob.matrix = prob.matrix
	} else if ( all(prob.matrix == last.prob.matrix) ) {
		cat('bootstrap PPM converged!\n')
		break
	} else {
		# mix the match-generated motif with the search motif
		# how...?
		# lets just try weighted average probabilities...
		prob.matrix = (1-opt$mixture) * prob.matrix + opt$mixture * last.prob.matrix
	}

	mname = sprintf('BS.%s.%s.%04d',opt$motf,opt$seqf,iter)
	# write new bootstrap motif to file, updating 'mfile' in the process so that it is used next
	mfile = write.MEME.motif(prob.matrix,bg.probs,mname)
}

if(opt$images) dev.off()
