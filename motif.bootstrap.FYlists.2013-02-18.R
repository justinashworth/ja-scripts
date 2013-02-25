#!/usr/bin/env Rscript

# search target sequence (e.g. genome) with a starting MEME motif,
#1 find matches (to do: optionally filtered to certain regions, e.g. promoters or ChIP overlaps?),
#2 build new motif based on actual matched sequences
#3 re-search with new [mixture] motif (bootstrap)
# repeat 1-3 until satisfied

# example
#./motif.bootstrap.R -m T.vol.FL3.MEME -s TSS-TTS.txt -n 100 -x 500 -b halo.bg.file -p 10 -i -X 0.25

source('motif.bootstrap.lib.R')

# INIT
library('getopt')
opt = getopt(matrix(c(
	'seqf'         , 's', 1, 'character',
	'niter'        , 'n', 1, 'integer',
	'maxhits'      , 'h', 1, 'integer',
	'mixture_start', 'X', 1, 'double',
	'mixture_end'  , 'x', 1, 'double',
	'finish_iters' , 'f', 1, 'integer',
	'bgfile'       , 'b', 1, 'character',
	'regions'      , 'r', 1, 'character',
	'genes'        , 'g', 1, 'character',
	'pseudocounts' , 'p', 1, 'double',
	'mt_start'     , 'M', 1, 'double',
	'mt_end'       , 'm', 1, 'double',
	'images'       , 'i', 0, 'logical',
	'help'         , 'H', 0, 'logical'
),ncol=4,byrow=T))

if( is.null(opt$seqf)         ) opt$seqf = 'halo.genome.fa'
if( is.null(opt$niter)        ) opt$niter = 50
#if( is.null(opt$niter)        ) opt$niter = 1
if( is.null(opt$maxhits)      ) opt$maxhits = 200
if( is.null(opt$mixture_start)) opt$mixture_start = 0.75
if( is.null(opt$mixture_end)  ) opt$mixture_end = 0.25
if( is.null(opt$finish_iters) ) opt$finish_iters = 10
if( is.null(opt$bgfile)       ) opt$bgfile = 'halo.bg.file'
#if( is.null(opt$regions)      ) opt$regions = 'tssflat.tsv'
# because we are generally comparing to ChIP data in direct upstream regions in other analyses(?),
# and because the TSS table is still considerably incomplete,
# we just look in regions directly upstream of the annotated ORF starts,
# rather than considering operonic promoters or true TSSs
if( is.null(opt$regions)      ) opt$regions = 'halo.gene.coords.tsv'
if( is.null(opt$pseudocounts) ) opt$pseudocounts = 1
if( is.null(opt$mt_start)     ) opt$mt_start = 0.0005
if( is.null(opt$mt_end)       ) opt$mt_end = 0.0001
if( is.null(opt$images)       ) opt$images = TRUE

print(opt)

# read/traverse Fang Yin's gene lists,
#		matching pred motifs to factors,
#			yeilding refit predicted motifs for each factor (under each condition?)

load('FY.Lrp.Lists.RData')

# TF translation table
tfcodes = names(neg.correlated.ChIP.targets)
lookup = read.delim('~/baliga/halo/Trh/sequence/VNG.codes')
tfnames = sapply(tfcodes, function(x){ lookup$name[ lookup$id==x ] })

# append SELEX motifs to lists
# TvFL3 <-> trh3, SsLrpB <-> trh6
tfcodes = c(as.character(tfcodes), 'VNG1816G', 'VNG1351G')
tfnames = c(as.character(tfnames), 'TvFL3', 'SsLrpB')

tftable = data.frame(code=tfcodes,name=tfnames)
print(tftable)

motfiles = list(
	'asnC'     = 'Halo.AsnC.MEME',
	'trh2'     = 'Halo.Trh2.MEME',
	'trh3'     = 'Halo.Trh3.MEME',
	'trh4'     = 'Halo.Trh4.MEME',
	'trh6'     = 'Halo.Trh6.MEME',
	'trh7'     = 'Halo.Trh7.MEME',
	'VNG1179C' = 'Halo.VNG1179C.MEME',
	'VNG1237C' = 'Halo.VNG1237C.MEME',
	'TvFL3'    = 'Tv.FL3.SELEX.MEME',
	'SsLrpB'   = 'Ss.LrpB.SELEX.MEME'
)

# make 'genomewide' list structure
regions = read.delim(opt$regions,as.is=T)
allgenes = unique(regions$canonical_Name)
genomewide = lapply(Lrp.ChIP.targets, function(x){allgenes})
names(genomewide) = names(Lrp.ChIP.targets)

inputs = list(
#	'genomewide'   = list(genomewide),
	'ChIP'         = list(Lrp.ChIP.targets),
#	'pos.ChIP'     = list(pos.correlated.ChIP.targets),
#	'neg.ChIP'     = list(neg.correlated.ChIP.targets),
	'pos.neg.ChIP' = list(pos.correlated.ChIP.targets, neg.correlated.ChIP.targets)
#	'pos.neg.cor'  = list(pos.correlated.genes, neg.correlated.genes)
#	'pos.ChIP.KO'  = pos.correlated.ChIP.targets.overlapwith.KO,
#	'neg.ChIP.KO'  = neg.correlated.ChIP.targets.overlapwith.KO,
#	'pos.genes.KO' = pos.correlated.genes.overlapwith.KO,
#	'neg.genes.KO' = neg.correlated.genes.overlapwith.KO
#	'pos.neg.ChIP.KO'  = c('pos.ChIP.KO','neg.ChIP.KO'),
#	'pos.neg.genes.KO' = c('pos.genes.KO','neg.genes.KO')
)

genelists = list()
for(input in names(inputs)){
	cat(input,'\n')
	genelists[[input]]=list()

	# pair TF names with MEME motifs
	for(tf in 1:nrow(tftable)){
		tfcode = as.character(tftable$code[tf])
		tfname = as.character(tftable$name[tf])
		motfile = as.character(motfiles[[tfname]])
		pair = paste(tfcode,motfile,sep='.')
		cat('\t',pair,'\n')

		# make gene list to search with motif
		genelists[[input]][[pair]] =
			unique( unlist( sapply(inputs[[input]], function(x){ x[[tfcode]] } ) ))
	}
}

cat('adapting predicted motifs\n')

results = list()

for(tf in 1:nrow(tftable)){
	tfcode = as.character(tftable$code[tf])
	tfname = as.character(tftable$name[tf])
	motfile = as.character(motfiles[[tfname]])
	pair = paste(tfcode,motfile,sep='.')
	cat('\t',pair,'\n')
	results[[pair]] = list()

	for(genelist in names(genelists)){
		cat(genelist,'\n')
		genes = genelists[[genelist]][[pair]]
		genes = unique( unlist(genes) )
		genes = genes[ !is.na(genes) ]
		if(length(genes)<1){
			cat('no genes!\n')
			next
		}
		cat(length(genes),'genes\n')
		print( head(genes) )

		prefix=paste(pair,'/',genelist,'/',sep='')
		cat(prefix,'\n')
		dir.create(prefix,recursive=T)

		result = motif_bootstrap(
			motfile=motfile,
			seqfile=opt$seqf,
			bgfile=opt$bgfile,
			regionsfile=opt$regions,
			genelist=genes,
			niter=opt$niter,
			maxhits=opt$maxhits,
			mixture_start=opt$mixture_start,
			mixture_end=opt$mixture_end,
			finish_iters=opt$finish_iters,
			pseudocounts=opt$pseudocounts,
			mt_start=opt$mt_start,
			mt_end=opt$mt_end,
			images=opt$images,
			prefix=prefix
		)
		save(result,file=sprintf('%sresult.RData',prefix))
		results[[pair]][[genelist]] = result

		# test for enrichment of matched genes in 'correlated' genes
		cor.genes = unique( unlist( genelists[['pos.neg.cor']][[pair]] ) )
		q = length( which(cor.genes %in% result$all.genes.matched) )
		m = length( result$all.genes.matched )
		n = length( result$all.genes ) - length( result$all.genes.matched )
		k = length( cor.genes )
		pv = phyper(q,m,n,k)
		results[[pair]][[genelist]][['fit.vs.cor.p']] = pv

		cat('fit', pair, genelist, 'vs. unfit cor genes:', q, m, n, k, pv, '\n')
	}

}

save(results,file='motif.bootstrap.results.RData')
