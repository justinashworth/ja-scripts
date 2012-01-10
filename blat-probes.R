jumps = function(x)
{
	jmps=c()
	x=sort(x)
	for ( i in 1:length(x)-1 ) {
		jmps[i]=x[i+1]-x[i]
	}
	return(jmps)
}

probeGaps = function(data)
{
	l=levels(factor(data$Tname))
	gaps=c()
	strands=c('+','-')
	for (i in 1:length(l)){
		level=l[i]
		for (j in 1:length(strands)){
			strand=strands[j]
			x=data$Tstart[which(data$Tname==level & data$strand==strand)]
			gaps=c(gaps,jumps(x))
		}
	}
	return(gaps)
}

analMultiMatch = function(data,maxmatch=1)
{
	#blat data in R-readable form
	f=factor(data$probe)
	t=table(f)
	print( paste('probes:',length(t)) )
	# number of probes matching multiple locations
	print( paste('multimatched:',length(t[which(t>maxmatch)]) ) )
	# plot matches per probe, ditributed by row order
	plot(t)
	# filter out mismatches
	return( data[ data$probe %in% t[which(t<=maxmatch)], ] )
}
