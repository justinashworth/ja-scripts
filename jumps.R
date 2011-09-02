jumps = function(x) {
	jmps=c()
	x=sort(x)
	for ( i in 1:length(x)-1 ) {
		jmps[i]=x[i+1]-x[i]
	}
	return(jmps)
}
