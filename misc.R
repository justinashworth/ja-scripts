# miscellaneous R hacks that may or may not be used or useful

pairwise =
	# faster than combn?
	# faster than for(){for(){}}?
	function(n)
{
	do.call('cbind',
		sapply(1:(n-1),function(i){
			sapply((i+1):n,function(j){c(i,j)})
		}))
}

peakseq =
	function(pos,seq,window=c(-200,200))
{
	start=max(0,pos+window[1])
	end=min(pos+window[2],nchar(seq))
	return(substring(seq,start,end))
}

peakseqs =
	# peaks: list of named vectors of peak positions (each list of peaks corresponds to a differently-named sequence region)
	# e.g.: list( 'Chr1'=vector([peak positions]), 'Chr2'=vector([peak positions]) )
	# seqs: list of named sequences
	# e.g.: list( 'Chr1'='ACGT...', 'Chr2'='ACGT...' )
	function(peaks,seqs,window=c(-200,200))
{
	pseqs=vector('list',length(peaks))
	names(pseqs)=names(peaks)
	for(seqname in names(peaks)){
		cat(seqname,'\n')
		if(! seqname %in% names(seqs)){
			cat('peaks for ',seqname,' have no corresponding sequence in provided list of sequences\n')
			next
		}
		pseqs[[seqname]] = sapply(peaks[[seqname]],peakseq,seqs[[seqname]],window)
		# 'unlist' will take care of the following later
#		names(pseqs[[seqname]]) = paste(seqname,peaks[[seqname]],sep='.')
		names(pseqs[[seqname]]) = peaks[[seqname]]
	}
	return(pseqs)
}

MeDiChI.peaks =
	# fits: a MeDiChI fit object
	function(fits,pval=0.05)
{
	# copy MeDiChI peaks into a simple list
	seqs=names(fits$fits.fin)
	peaks=list()
	for(seq in seqs){
		seqname=seq
		if(seq=='HALCHR'){seqname='Chr'}
		peaks[[seqname]]=fits$fits.fin[[seq]][[1]]$coeffs.w.p.values
	}
	# filter by pval
	peaks=sapply(peaks,function(x){x[x[,3]<pval,1]})
	# round peak positions
	peaks=sapply(peaks,round)
	return(peaks)
}

load.sequence =
	# convenience function for loading fasta sequence into simple list format
	function(file='halo.genome.sequence.fa')
{
	library(Biostrings)
	seqs=readFASTA(file)
	# convert to better list structure
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

