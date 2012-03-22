### Justin Ashworth 2010
### Institute for Systems Biology

### Data processing and plotting for manual and automatic measurements of Thalassiosira bioreactor culture

dec.hours =
	function(x)
	# H:M:S to decimal hours
{
	s=strsplit(as.character(x),':')
	as.numeric(s[[1]][1]) + as.numeric(s[[1]][2])/60 + as.numeric(s[[1]][3])/60/60
}

convert.manual.time =
	function(manual.data)
{
	manual.data$date.time=strptime(manual.data$date.time,'%m/%d/%y %I:%M %p')
	manual.data$expt.time=sapply(manual.data$expt.time,dec.hours)
	return(manual.data)
}

unfactor =
	function(x){levels(x)[as.numeric(x)]}

delim.expand =
	function(x,y,delim=';')
	# expand paired data in which y has multiple delimited values for individual x values
	# works on (and filters/extends) a pair of corresponding vectors
	# simplifies bookkeeping and data entry
{
	if(is.factor(x)){x=unfactor(x)}
	if(is.factor(y)){y=unfactor(y)}
	for(i in grep(delim,y)){
		sp=unlist(strsplit(as.character(y[i]),delim))
		for(s in sp){
			#cat(x[i],s,'\n')
			x=c(x,x[i])
			y=c(y,s)
		}
	}
	xy=cbind(x,y)
	return(xy[!grepl(delim,y) & !y=="",])
}

delim.mean.sd =
	function(x,delim=';')
	# not vectorized: use an apply method
{
	sp=unlist(strsplit(as.character(x),delim))
	return(c(mean(as.numeric(sp),na.rm=T),sd(as.numeric(sp),na.rm=T)))
}

growthrate.specific =
	function(x,y)
{
	l=length(x)
	if(l != length(y)){
		cat('incompatible x and y vectors\n')
		return()
	}
	return( log(y[2:l] / y[1:l-1]) / (x[2:l] - x[1:l-1]) )
}

growthrate =
	function(x,y)
{
	l=length(x)
	if(l != length(y)){
		cat('incompatible x and y vectors\n')
		return()
	}
	return((y[2:l]-y[1:l-1])/(x[2:l]-x[1:l-1]))
}

#doublingtime =
#	log(2)/log(1+growthrate)


