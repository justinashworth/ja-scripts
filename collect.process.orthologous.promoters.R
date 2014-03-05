# Justin Ashworth 2014
# Institute for Systems Biology

# Original code for multi-species orthologous promoter analysis for conserved non-coding (cis-regulatory) motifs
# 1. protein BLAST to find orthologous proteins between query genome and all subjects
# 2. extraction of the promoter region sequences of these orthologous genes in multiple species
# 3. MEME for conserved motifs
# 4. [optional] MAST back to the promoter regions of single query/reference genome
# [5. refactor, optimize for speed, including multiprocessing... and/or port to c++ using NCBI toolbox libraries which could be 10 or 100 times faster?]
# [6. combine with cMonkey]

################################################################################
# FUNCTIONS AND GLOBALS (most/all global data should be UPPERCASE)

# download genome data (gff, faa, fna)
# This code is driven by fetching corresponding data from faa, fna, and gff files from NCBI.
# There may be some alternative scheme involving the elegant processing of ASN.1 or Genbank files, but the basic idea is the same
NCBI_downloads =
	function(tab='halobacteria.tax')
{
	d=read.delim(tab)
	suff = c('faa','fna','gff')

	apply(d, 1, function(x){
		cat(x[['taxid']],'\n');
		sapply(suff, function(sf){
			cat(sf,' ')
			cmd = sprintf('wget "%s*%s"', x[['NCBIftpdir']], sf)
			system(cmd)
		})
		cat('\n')
	})
}

# preload all genome sequences and gffs as globals, for speed/efficiency
require(Biostrings)
system('ls *.fna > allfnas')
cat('Preloading all genome sequences...')
fnas = readLines('allfnas')
ALLGENOMESEQS = DNAStringSet( sapply(fnas, function(x){ readDNAStringSet(x)[[1]] }) )
names(ALLGENOMESEQS) = gsub('.fna','',fnas)
cat('done\n')
# GFFs
system('ls *.gff > allgffs')
cat('Preloading all GFFs...')
gffs = readLines('allgffs')
colnames = c('genome','source','feature','start','end','score','strand','frame','attr')
ALLGFFS = lapply(gffs, function(x){ d = read.delim(x,header=F,comment.char='#'); colnames(d)=colnames; d })
names(ALLGFFS) = gsub('.gff','',gffs)
cat('done\n')

# use of this parameter is what excludes both i) coding sequences 5' of the promoters and ii) genes internal to operons
MINPROMLEN = 20

# produce promoter sequences for all hits
prom_seqs =
	# note that using a positive value for reg[2] includes protein-coding sequence. This is nice for illustrating that promoters are properly aligned in alignments. It's very bad for MEME though, so there is some extra logic in this function and in run_meme() to trim this out for meme and filtering purposes.
	function(protids,reg=c(-100,20))
{
	require(Biostrings)
	seqs = sapply(protids, function(x){

		# find protein in a gff file
		# lo-tech grep for protein id in gffs (replace with something fancy?)
		gff = system( sprintf('grep %s *gff',x), intern=T)
		# parse gff coords
		gff = unlist(strsplit(gff,"\t"))
		refseq = gff[1]
		refseq = as.character(unlist(strsplit(refseq,'\\.'))[1])
		gstart = as.integer(gff[4])
		gend = as.integer(gff[5])
		strand = gff[7]
		#cat(refseq,gstart,gend,strand,'\n')

		# define sequence region
		regstart = gstart + reg[1]
		regend = gstart + reg[2]
		if(strand=='-'){
			regstart = gend - reg[2]
			regend = gend - reg[1]
		}
		#cat(regstart,regend,'\n')

		## start code to avoid coding sequences 5' of the promoter, as well as genes internal to operons ##
		# find nearest 5' CDS, trim 'promoter' region
		if(!refseq %in% names(ALLGFFS)){
			cat('PROBLEM:',refseq,'not in GFFs!!!\n')
			return(NA)
		}
		gff = ALLGFFS[[refseq]]

		if(strand=='+'){
			# not assuming GFF is ordered (though it should be)
			prevcds = gff[gff$end <= gstart,]
			if(nrow(prevcds)>0){
				prevcds = prevcds[ order(prevcds$end,decreasing=T)[1], ]
				#cat('previous cds is', as.character(prevcds),'\n')
				if(prevcds[['end']] > regstart){
					#cat('set regstart',regstart)
					regstart=prevcds[['end']]
					newlen = regend-regstart
					#cat('to',regstart,'for new promoter length',newlen,'\n')
					# skip if too short/nonexistent (e.g. inside operon)
					# note that reg[2] must be added to MINPROMLEN since we are adding/subtracting reg[2] from the length of promoter sequences
					if(newlen < (MINPROMLEN+reg[2])){
						cat('Skipping promoter trimmed to length',regend-regstart-reg[2],'\n')
						return(NA)
					}
				}
			}

		}else if(strand=='-'){
			nextcds = gff[gff$start >= gend,]
			if(nrow(nextcds)>0){
				nextcds = nextcds[ order(nextcds$start)[1], ]
				#cat('next cds is', as.character(nextcds),'\n')
				if(nextcds[['start']] < regend){
					#cat('set regend',regend)
					regend=nextcds[['start']]
					newlen = regend-regstart
					#cat('to',regend,'for new promoter length',newlen,'\n')
					# skip if too short/nonexistent (e.g. inside operon)
					# note that reg[2] must be added to MINPROMLEN since we are adding/subtracting reg[2] from the length of promoter sequences
					if(newlen < (MINPROMLEN+reg[2])){
						cat('Skipping promoter trimmed to length',regend-regstart-reg[2],'\n')
						return(NA)
					}
				}
			}
		}

		# extract promoter region sequence from genome
		if(regstart<1) {
			cat('Warning!!! regstart<1 for',x,refseq,gstart,gend,strand,regstart,regend,'\n')
			regstart = 1
		}
		if(!refseq %in% names(ALLGENOMESEQS)){
			cat('PROBLEM:',refseq,'not in genome sequences!!!\n')
			return(NA)
		}
		seqlen = nchar(ALLGENOMESEQS[[refseq]])
		if(regend>seqlen) {
			cat('Warning!!! regend is out of bounds for',x,refseq,seqlen,gstart,gend,strand,regstart,regend,'\n')
			regend = seqlen
		}

		seq = ALLGENOMESEQS[[refseq]][regstart:regend]
		if(strand=='-'){
			seq = reverseComplement(seq)
		}
		seq
	})
	names(seqs) = protids
	return(seqs)
}

MIN_MEME_SEQS = 5

run_meme =
	function(promseqs,name,prom_reg=c(0,0),bfile=NULL)
{
	meme_cmd = 'meme_4.9.0 %s -dna -oc %s -nmotifs 2 -minw 10 -maxw 16 -nostatus'
	cat('MEMEing for',name,'\n')

	if(length(promseqs)<2){
		cat('Less than 2 sequences for',name,'(skipping)\n')
		return()
	}

	seqf = sprintf('promseqs/%s.fa',name)
	promseqs = DNAStringSet(promseqs)
	writeXStringSet(promseqs, seqf)
	# do a MUSCLE alignment while we're at it
	# FASTA (better support for various programs, e.g. re-input to muscle itself)
	muscle = sprintf('muscle3.8.31_i86darwin64 -in %s -out promseqs/`basename %s`.aln.fa -quiet', seqf, seqf)
	system(muscle)
	# ClustalW format, more human-readable
	muscle = sprintf('muscle3.8.31_i86darwin64 -in %s -out promseqs/`basename %s`.aln.clw -quiet -clw', seqf, seqf)
	system(muscle)

	# remove duplicate sequences (likely an artifact of BLAST redundancies, duplicate gene models, etc)
	nseqs = length(promseqs)
	promseqs = unique(promseqs)
	if(length(promseqs)<nseqs){
		cat('removed',nseqs-length(promseqs),'duplicate promoter sequences\n')
	}

	if(length(promseqs)<MIN_MEME_SEQS){
		cat('not enough unique orthologous sequences for',name,'to compute a motif\n')
		return()
	}

	# trim out coding sequence for MEME if it is indicated to exist according to prom_reg
	if(prom_reg[2]>0){
		promseqs = DNAStringSet( sapply(promseqs, function(x){ x[ 1:(nchar(x)-prom_reg[2]) ] }) )
		writeXStringSet(promseqs, seqf)
	}
	cmd = sprintf(meme_cmd, seqf, sprintf('MEME/%s',name))
	if(!is.NULL(bfile) & file.exists(bfile)) cmd = sprintf('%s -bfile %s',cmd,bfile)
	system(cmd)
}

run_mast =
	function(name,targets='HalNRC-1.proms.fa',bfile='halo.bg.file',ev='100',mt='0.0001')
{
	memefile = sprintf('MEME/%s/meme.txt',name)
	if(!file.exists(memefile)) return()
	mast_cmd = 'mast_4.9.0 %s %s -bfile %s -oc %s -ev %s -mt %s -nostatus'
	cat('MASTing for',name,'vs',targets,'\n')
	cmd = sprintf(mast_cmd, memefile, targets, bfile, sprintf('MAST/%s_%s',name,targets), ev, mt)
	system(cmd)
}

make_BLAST_db =
	function(subjectdb_name,overwrite=F)
{
	if( !file.exists(sprintf('%s.psq',subjectdb_name)) | overwrite ){
		# makeblastdb for all proteins (subject db)
		system( sprintf('cat *faa > %s.faa', subjectdb_name) )
		system( sprintf('makeblastdb -in %s.faa -dbtype prot -title %s -out %s', subjectdb_name, subjectdb_name, subjectdb_name ) )
	}
}

# END FUNCTIONS AND GLOBALS
################################################################################
################################################################################
# PREPARE/PERFORM BLAST

# UNCOMMENT THIS only do this once!!
#NCBI_downloads()

# combine proteins from query genome
#NC_002607 Halobacterium sp. NRC-1 chromosome
#NC_001869 Halobacterium sp. NRC-1 plasmid pNRC100
#NC_002608 Halobacterium sp. NRC-1 plasmid pNRC200
query_accs = c('NC_001869','NC_002607','NC_002608')
prots = 'NRC-1'
system( sprintf('cat %s.faa > %s.faa', paste(query_accs,collapse='.faa '), prots) )
BLASTdb = 'halobacteria'

BLAST_result = sprintf('%s.BLAST',prots)
overwrite = F
# only do this once (time-consuming)
if(!file.exists(BLAST_result) | overwrite) {
	# make BLAST db for all genomes' proteins in existing directory (*faa)
	make_BLAST_db(BLASTdb)
	# run BLASTp for all query proteins vs all genomes' proteins
	blastcmd = sprintf('blastp -db %s -query %s.faa -out %s -outfmt 7 -evalue 1e-10', BLASTdb, prots, BLAST_result)
	cat(blastcmd,'...\n')
	system(blastcmd)
} else {
	cat('Using old precomputed BLAST results in current directory\n')
}

# extract promoter sequences for orthologous genes
cat('Loading BLAST results table...')
BLAST = read.delim(BLAST_result,comment.char='#',header=F)
# BLAST fields:
# query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
names(BLAST) = c('qry', 'sbj', 'ident', 'len', 'mis', 'gap', 'qstart', 'qend', 'sstart', 'send', 'eval', 'bits')
cat('done\n')

# parse/clean query/subject ids down to protein accession numbers
BLAST$qry = do.call(rbind, strsplit(as.character(BLAST$qry),'\\|'))[,4]
BLAST$sbj = do.call(rbind, strsplit(as.character(BLAST$sbj),'\\|'))[,4]

# END PREPARE/PERFORM BLAST
################################################################################
################################################################################
# The sequence collection routine

# TEST etfB
#BLAST = BLAST[ BLAST$qry=='NP_280807.1', ]

cat('Loading query proteins...')
QUERY_PSEQS = readDNAStringSet( sprintf('%s.faa',prots) )
names(QUERY_PSEQS) = do.call(rbind, strsplit(names(QUERY_PSEQS),'\\|'))[,4]
cat('done\n')

collect_orthologous_promoters =
	# partially subroutined...
	function(queries, min_ident, min_length, min_coverage, prom_reg)
{
	cat('Collecting othologous promoter sequences...\n')
	all_promseqs = lapply(queries, function(qry){
		cat(qry,'\n')
		qry_length = length(QUERY_PSEQS[[qry]])
		hits = BLAST[ BLAST$qry==qry & BLAST$ident>=min_ident & BLAST$len>=min_length & (BLAST$len/qry_length)>=min_coverage, ]
		cat(nrow(hits),'BLASTp hits\n')
		promseqs = prom_seqs(hits$sbj,reg=prom_reg)
		nNAs = length(which(is.na(promseqs)))
		if(nNAs>0) cat(nNAs,'promseqs were NAs\n')
		promseqs = promseqs[ !is.na(promseqs) ]
		promseqs
	})
	names(all_promseqs) = queries
}

################################################################################
################################################################################
# most user params and function calls below

# choose high identity hits (>=50%) with length >= 50 amino acids and over 70% query coverage
MIN_IDENT = 50
MIN_LENGTH = 50
MIN_COVERAGE = 0.7
PROM_REG=c(-100,20)

param_tag = sprintf('mid%i_mln%i_mcv%g_pr%i%i',MIN_IDENT,MIN_LENGTH,MIN_COVERAGE,PROM_REG[1],PROM_REG[2])

preserve = F
if(!preserve){
	system('find MEME -delete')
	system('find MAST -delete')
	system('find promseqs -delete')
}

system('mkdir MEME')
system('mkdir MAST')
system('mkdir promseqs')

# uncomment below to run actual procedure
query_prots = unique(BLAST$qry)
all_promseqs = collect_orthologous_promoters(query_prots,min_ident=MIN_IDENT,min_length=MIN_LENGTH,min_coverage=MIN_COVERAGE,prom_reg=PROM_REG)
save(all_promseqs,file=sprintf('OrthologousPromoterSeqs.%s.RData',param_tag))

# or run MEME on previously collected orthologous promoters
#load('OrthologousPromoterSeqs.RData')

# add protein names and species-specific gene ids to make browsing easier
ptable = read.delim('halo.proteins.tsv',as.is=T)
pnames = ptable$name
names(pnames) = ptable$acc
ssgids = ptable$vng
names(ssgids) = ptable$acc
names(all_promseqs) = paste(names(promseqs),pnames[names(promseqs)],ssgids[names(promseqs)],sep='_')

MAST = F
cat('Running MEME/MAST...\n')
for(prot in names(all_promseqs)){
	#cat(prot,'\n')
	# MEME the orthologous promoter sequences
	# feed 'prom_reg' to this function to avoid finding protein-coding motifs
	run_meme(all_promseqs[[prot]],prot,prom_reg=PROM_REG)
	# MAST the orthologous motifs against the query/'reference' genome
	if(MAST) run_mast(prot)
}
cat('Done\n')
