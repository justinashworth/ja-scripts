#!/usr/bin/env Rscript
# MAST target sequence (e.g. genome) with a starting MEME motif,
#1 find matches (optionally filtered to certain regions, e.g. promoters),
#2 build new motif based on actual matched sequences (mixture?)
#3 re-MAST with new motif (bootstrap)
# repeat 1-3 until...?

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
	'bfile'        , 'b', 1, 'character',
	'pseudocounts' , 'p', 1, 'integer',
	'png'          , 'g', 0, 'logical'
),ncol=4,byrow=T))

if( is.null(opt$motf) | is.null(opt$seqf) ) {
	cat('args: -m <motif_file> -s <fasta_seq>\n')
	opt$motf = 'T.vol.FL3.MEME'
	opt$seqf = 'halo.genome.fa'
	opt$bfile = 'halo.bg.file'
#	q('no')
}

if( is.null(opt$niter) ) opt$niter = 10
if( is.null(opt$maxhits) ) opt$maxhits = 500
if( is.null(opt$pseudocounts) ) opt$pseudocounts = 1
if( is.null(opt$png) ) opt$png = FALSE

bg.probs = list()

if( is.null(opt$bfile) ) {
	bg.probs = list(
		'A'=0.25,
		'C'=0.25,
		'G'=0.25,
		'T'=0.25
	)
} else {
	for( line in readLines(opt$bfile) ) {
		fields = unlist(strsplit(line,'[ \t]'))
		bg.probs[fields[1]] = as.numeric(fields[2])
	}
}

mast.args = '-hit_list'
if( !is.null(opt$bfile) ) { mast.args = paste(mast.args,'-bfile',opt$bfile) }

cat('background probabilities:\n')
print(paste(unlist(bg.probs)),sep='')

cat('NOTE: currently stripping letter info in bg probs and assuming implicit A,C,G,T ordering!\n')

if(opt$png) png('BS.motif.%d.png')

mfile = opt$motf
last.prob.matrix = NULL
for( iter in 1:opt$niter ) {
	cat('iter ',iter,'\n',sep='')
	cat('MASTing motif from file "',mfile,'" in sequence from fasta file "',opt$seqf,'"\n',sep='')
	hitfile = paste('hitfile.tmp.',iter,sep='')
	mast.cmd = paste('mast',mfile,opt$seqf,mast.args,'>',hitfile)
	cat(mast.cmd,'\n')
	system(mast.cmd)
	hits = read.table(hitfile)

	# sort by pval and cull
	hits=hits[ order(hits[,6])[ 1:min(nrow(hits),opt$maxhits) ] , ]

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
		cat('bootstrap ppm converged!\n')
		break
	}
	mname = paste('BS',opt$motf,opt$seqf,iter,sep='.')
	mfile = write.MEME.motif(prob.matrix,bg.probs,mname)
}

if(opt$png) dev.off()
