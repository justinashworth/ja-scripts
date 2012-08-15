#!/usr/bin/env Rscript

# MAST target sequence (e.g. genome) with a starting MEME motif,
#1 find matches (to do: optionally filtered to certain regions, e.g. promoters or ChIP overlaps?),
#2 build new motif based on actual matched sequences
#3 re-MAST with new [mixture] motif (bootstrap)
# repeat 1-3 until satisfied

# example
#./MAST.bootstrap.R -m T.vol.FL3.MEME -s TSS-TTS.txt -n 100 -x 500 -b halo.bg.file -p 10 -i -X 0.25

source('MAST.bootstrap.lib.R')

# INIT
library('getopt')
opt = getopt(matrix(c(
	'motf'         , 'm', 1, 'character',
	'seqf'         , 's', 1, 'character',
	'niter'        , 'n', 1, 'integer',
	'maxhits'      , 'x', 1, 'integer',
	'mixture'      , 'X', 1, 'double',
	'bgfile'       , 'b', 1, 'character',
	'regions'      , 'r', 1, 'character',
	'genes'        , 'g', 1, 'character',
	'pseudocounts' , 'p', 1, 'integer',
	'images'       , 'i', 0, 'logical',
	'help'         , 'h', 0, 'logical'
),ncol=4,byrow=T))

if( is.null(opt$niter)        ) opt$niter = 10
if( is.null(opt$maxhits)      ) opt$maxhits = 500
if( is.null(opt$mixture)      ) opt$mixture = 0.25
if( is.null(opt$pseudocounts) ) opt$pseudocounts = 1
if( is.null(opt$images)       ) opt$images = FALSE

if( is.null(opt$motf) | is.null(opt$seqf) | !is.null(opt$help) ) {
	cat('args: -m <motif_file> -s <fasta_seq> [-[n|-niter] <%d>] [-[x|-maxhits] <%d>] [-[b|-bgfile] <bgfile>] [-[p|-pseudocounts] <%d>] [-[g|-images]]\n')
	opt$motf = 'T.vol.FL3.MEME'
	opt$seqf = 'halo.genome.fa'
	opt$bgfile = 'halo.bg.file'
	#opt$regions = 'halo.gene.coords.tsv'
	opt$regions = 'TSS-TTS.txt'
#	q('no')
}

if( opt$mixture < 0 | opt$mixture > 1 ){
	cat('invalid mixture ratio\n')
	q('no')
}

print(opt)

motif_bootstrap(opt$motf, opt$seqf, opt$bgfile, opt$regions, opt$niter, opt$maxhits, opt$mixture, opt$pseudocounts, opt$images)
