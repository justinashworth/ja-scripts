probMatrixFromFasta =
	function(fastafile,filter='')
{
	require(Biostrings)
	seqs=read.DNAStringSet(fastafile)
	if(filter != '') seqs = gsub(filter,'',seqs)
	ppm=consensusMatrix(seqs,as.prob=T,baseOnly=T)
	ppm=as.matrix(ppm)
	return(ppm)
}

seqlogoFromFasta =
	function(fastafile,plot=F)
{
	require(seqLogo)
	ppm = probMatrixFromFasta(fastafile)
	# limit to first four rows (ACGT) (drops fifth row, 'other')
	ppm=ppm[seq(1,4),]
	# hardcoded arbitrary subsequence for 2e1c DNA (otherwise rosetta 'design_mutations' DNA .fa files read through both strands)
	#ppm=ppm[,seq(2,16)]
	# normalize (if necessary)
	ppm=t(t(ppm)*1/colSums(ppm))
	if(plot){seqLogo(ppm)}
#	if(plot){seqLogo(ppm,ic.scale=F)}
	return(ppm)
}

dis.ED = function(P,Q)
	# Euclidian distance
{
	sqrt( sum((P-Q)**2) )
}

dis.KL = function(P,Q)
	# Kullback-Leibler divergence (column vs. column)
{
	sum(P * log(P/Q))
}

dis.JS = function(P,Q)
	# Jensen-Shannon divergence (column vs. column)
{
	M = 0.5*(P+Q)
	0.5 * ( dis.KL(P,M) + dis.KL(Q,M) )
}

ppm.distance.ED = function(ppm1,ppm2)
{
	motif_distance = 0
	for(col in 1:ncol(ppm1)){
		motif_distance = motif_distance + dis.ED( ppm1[,col], ppm2[,col] )
	}
	return(motif_distance)
}

ppm.distance.KL.symm =
	# works on two position probability matrices
	# bases in rows, positions in columns. Note correspondence check should have happened upstream
	# computes symmetric Kullback-Leibler distance between them
	function(ppm1,ppm2)
{
	motif_distance = 0
	for(col in 1:ncol(ppm1)){
		motif_distance = motif_distance + 0.5 * ( dis.KL(ppm1[,col], ppm2[,col]) + dis.KL(ppm2[,col], ppm1[,col]) )
	}
	return(motif_distance)
}

ppm.distance.ALLR =
	# works on two aligned position probability matrices
	# bases in rows, positions in columns. Note correspondence check should have happened upstream
	# computes Wang and Stormo average log likelihood ratio (ALLR) statistic
	# because the input is probabilities and not counts, total column counts are assumed to be constant between columns and motifs, and thus base probabilities are used in place of base counts
	function(ppm1,ppm2,bg=rep(0.25,4),nonzero=1e-2)
{
	motif_dist = 0
#	maxdis = 10
	zeroflag=0
	for(col in 1:ncol(ppm1)){
		ALLR = 0
		for(base in 1:nrow(ppm1)){
			p1 = ppm1[base,col]
			p2 = ppm2[base,col]
			ballr=0
			if( (p1==0 | p2==0) & zeroflag){
				cat('zeros encountered. Substituting',nonzero,'This can be avoided by using pseudocounts\n')
				zeroflag=0
			}
			p1 = max(p1,nonzero)
			p2 = max(p2,nonzero)
			ballr = (p1 * log(p1/bg[base]) + p2 * log(p2/bg[base])) / (p1+p2)
#			ballr = min(ballr,maxdis)
#			cat('column',col,'row',base,'ballr',ballr,'\n',sep=' ')
			ALLR = ALLR + ballr
		}
#		cat('column',col,'ALLR',ALLR,'\n',sep=' ')
		motif_dist = motif_dist + ALLR
	}
	return(motif_dist)
}

ppm.similarity.BLiC.Yanover =
	function(ppm1,ppm2,bg,...)
{
	motif_distance = 0
	for(col in 1:ncol(ppm1)){
		col_dis = dis.JS( ppm1[,col]+ppm2[,col], bg) - dis.JS( ppm1[,col], ppm2[,col] )
		motif_distance = motif_distance + col_dis
	}
	return(motif_distance)
}

# toying with Dirichlet distributions after reading Habib 2008
dirichlet.norm = function(alphas)
{
	gamma(sum(alphas)) / prod(gamma(alphas))
}

dirichlet =
	function(x,alphas)
{
	norm = dirichlet.norm(alphas)
	cat('norm',norm,'\n')
	if(length(x) != length(alphas)){
		cat('ERROR a,alphas unequal length\n')
		return(Inf)
	}
	val = 1
	for(i in 1:length(x)){
		val = val * (x[i]**(alphas[i]-1))
	}
	return(val * norm)
}

ppm.similarity.BLiC =
	# computes Bayesian Liklihood 2-Component (BLiC) score (Habib et al. 2008)
	# works on two position probability matrices
	# bases in rows, positions in columns. Note correspondence check should have happened upstream
	# uses dirichlet prior ('P12') to represent the "common source distribution"
	# in the simple case this is like adding pseudocounts??
	function(ppm1,ppm2,bg=rep(0.25,4),nonzero=1e-2,param=NULL)
{
	if(param==NULL){
		npos = ncol(ppm1)
		param=rep(1,pos)
	}

	# INCOMPLETE
	# WHAT IS P12??

	motif_distance = 0
#	maxdis = 10
	zeroflag = 0
	for(col in 1:ncol(ppm1)){
		column_distance = 0
		for(base in 1:nrow(ppm1)){
			p1 = ppm1[base,col]
			p2 = ppm2[base,col]
			if( (p1==0 | p2==0) & zeroflag){
				cat('zeros encountered. Substituting',nonzero,'This can be avoided by using pseudocounts\n')
				zeroflag=0
			}
			p1 = max(p1,nonzero)
			p2 = max(p2,nonzero)

			# "score is composed of two components: the first measures whether the two motifs were generated from a common distribution, while the second reflects the distance of that common distribution from the background"
			# BLiC = log( P(m1,m2|common-source)/P(m1,m2|independent-source) )
			#        + log( P(m1,m2|common-source)/P(m1,m2|background) )
			# P12 is a dirichlet mixture prior
			# not sure if this needs to fake 'counts' from input probability matrices or not
#			2 * (n1+n2) * log(P12[base]) - n1 * log(p1) - n2 * log(p2) - (n1+n2) * log(bg[base])
			2 * (p1+p2) * log(P12[base]) - p1 * log(p1) - p2 * log(p2) - (p1+p2) * log(bg[base])

#			cat('column',col,'row',base,'base_lr',base_lr,'\n',sep=' ')
			column_distance = column_distance + base_lr
		}
#		cat('column',col,'column_distance',column_distance,'\n',sep=' ')
		motif_distance = motif_distance + column_distance
	}
	return(motif_distance)
}

ppm.dis.matrix =
	function(matrixlist,disfunc='KL',bg=rep(0.25,4),...)
{
	if(0 %in% bg){
		cat('ERROR, no background probabilities may be zero\n')
		return(0)
	}
	names=names(matrixlist)
	l=length(matrixlist)
	dmat = matrix(0,nrow=l,ncol=l,dimnames=list(names,names))
	for(i in 1:l){
		for(j in 1:l){
			ppm1 = matrixlist[[i]]
			ppm2 = matrixlist[[j]]

			if(ncol(ppm1) != ncol(ppm2)){
				cat('ERROR: column mismatch\n')
				next
			}
			if(nrow(ppm1) != nrow(ppm2)){
				cat('ERROR: column mismatch\n')
				next
			}

			if(disfunc=='ED'){
				dmat[i,j] = ppm.distance.ED(ppm1,ppm2,...)
			} else if(disfunc=='ALLR'){
				dmat[i,j] = ppm.distance.ALLR(ppm1,ppm2,bg,...)
			} else if(disfunc=='KL'){
				dmat[i,j] = ppm.distance.KL.symm(ppm1,ppm2,...)
			} else if(disfunc=='BLiC.inv'){
				# BLiC-like similarity metric from Yanover 2011
				# sign inverted (for uniformity with distance metrics)
				dmat[i,j] = -1 * ppm.similarity.BLiC.Yanover(ppm1,ppm2,bg,...)
			} else {
				dmat[i,j] = ppm.distance.KL.symm(ppm1,ppm2,...)
			}
		}
	}
	return(dmat)
}

ppm.dis.test =
	function(disfunc='KL',bg=rep(0.25,4),...)
{
	A = c(1,0,0,0)
	C = c(0,1,0,0)
	G = c(0,0,1,0)
	T = c(0,0,0,1)
	ppms = list(
		 'AAAA' = matrix(c(A,A,A,A),4)
		,'AAAT' = matrix(c(A,A,A,T),4)
		,'AATT' = matrix(c(A,A,T,T),4)
		,'ATTT' = matrix(c(A,T,T,T),4)
		,'TTTT' = matrix(c(T,T,T,T),4)
	)
	return(ppm.dis.matrix(ppms,disfunc,bg,...))
}

self.distance.to.NA = function(mat){
	rns = rownames(mat)
	cns = colnames(mat)
	if(rns!=cns){return(mat)}
	for(rn in rns){
		mat[rns=rn,cns=rn]=NA
	}
	return(mat)
}
