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
	'seqf'         , 's', 1, 'character',
	'niter'        , 'n', 1, 'integer',
	'maxhits'      , 'x', 1, 'integer',
	'mixture'      , 'X', 1, 'double',
	'bgfile'       , 'b', 1, 'character',
	'regions'      , 'r', 1, 'character',
	'genes'        , 'g', 1, 'character',
	'pseudocounts' , 'p', 1, 'integer',
	'mt'           , 'm', 1, 'double',
	'images'       , 'i', 0, 'logical',
	'help'         , 'h', 0, 'logical'
),ncol=4,byrow=T))

if( is.null(opt$seqf)         ) opt$seqf = 'halo.genome.fa'
#if( is.null(opt$niter)        ) opt$niter = 10
if( is.null(opt$niter)        ) opt$niter = 1
if( is.null(opt$maxhits)      ) opt$maxhits = 250
if( is.null(opt$mixture)      ) opt$mixture = 0.7
if( is.null(opt$bgfile)       ) opt$bgfile = 'halo.bg.file'
if( is.null(opt$regions)      ) opt$regions = 'tssflat.tsv'
if( is.null(opt$pseudocounts) ) opt$pseudocounts = 1
if( is.null(opt$mt)           ) opt$mt = 0.0004
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

FY.listtypes = list(
	'pos.ChIP'     = pos.correlated.ChIP.targets,
	'neg.ChIP'     = neg.correlated.ChIP.targets,
	'pos.ChIP.KO'  = pos.correlated.ChIP.targets.overlapwith.KO,
	'neg.ChIP.KO'  = neg.correlated.ChIP.targets.overlapwith.KO,
	'pos.genes.KO' = pos.correlated.genes.overlapwith.KO,
	'neg.genes.KO' = neg.correlated.genes.overlapwith.KO
)

FY.combine = list(
	'pos.neg.ChIP'     = c('pos.ChIP','neg.ChIP'),
	'pos.neg.ChIP.KO'  = c('pos.ChIP.KO','neg.ChIP.KO'),
	'pos.neg.genes.KO' = c('pos.genes.KO','neg.genes.KO')
)

# create 'all' (pos+neg) gene lists
listtypes = list()

cat('creating \'all\' lists (pos+neg)\n')
for(listtype in names(FY.combine)){
	cat(listtype,'\n')
	listtypes[[listtype]]=list()
	for(tf in 1:nrow(tftable)){
		tfcode = as.character(tftable$code[tf])
		tfname = as.character(tftable$name[tf])
		motfile = as.character(motfiles[[tfname]])
		pair = paste(tfcode,motfile,sep='.')
		cat('\t',pair,'\n')
		listtypes[[listtype]][[pair]] =
			unique( unlist( sapply(
				FY.combine[[listtype]],
				function(x){ unlist(FY.listtypes[[x]][[tfcode]]) } )
			))
	}
}

cat('adapting predicted motifs\n')
for(tf in 1:nrow(tftable)){
	tfcode = as.character(tftable$code[tf])
	tfname = as.character(tftable$name[tf])
	motfile = as.character(motfiles[[tfname]])
	pair = paste(tfcode,motfile,sep='.')
	cat('\t',pair,'\n')

	for(listtype in names(listtypes)){
		cat(listtype,'\n')
		genelist = listtypes[[listtype]][[pair]]
		genelist = unique( unlist(genelist) )
		genelist = genelist[ !is.na(genelist) ]
		if(length(genelist)<1){
			cat('no genes!\n')
			next
		}
		cat(length(genelist),'genes\n')
		print( head(genelist) )

		prefix=paste(pair,'/',listtype,'/',sep='')
		dir.create(prefix,recursive=T)

		motif_bootstrap(
			motfile=motfile,
			seqfile=opt$seqf,
			bgfile=opt$bgfile,
			regionsfile=opt$regions,
			genelist=genelist,
			niter=opt$niter,
			maxhits=opt$maxhits,
			mixture=opt$mixture,
			pseudocounts=opt$pseudocounts,
			mt=opt$mt,
			images=opt$images,
			prefix=prefix
		)
	}

}

