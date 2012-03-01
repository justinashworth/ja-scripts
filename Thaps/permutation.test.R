# column definitions

dark = c(
	'X400.Day1.Dk',
	'X400.Day2.Dk',
	'X400.Day3.Dk',
	'X400.Day4.Dk',
	'X400.Day5.Dk',
	'X800.Day1.Dk',
	'X800.Day2.Dk',
	'X800.Day3.Dk',
	'X800.Day4.Dk',
	'X800.Day5.Dk'
)

light = c(
	'X400.Day1.Lt',
	'X400.Day2.Lt',
	'X400.Day3.Lt',
	'X400.Day4.Lt',
	'X400.Day5.Lt',
	'X800.Day1.Lt',
	'X800.Day2.Lt',
	'X800.Day3.Lt',
	'X800.Day4.Lt'
)

early = c(
	'X400.Day1.Dk',
	'X400.Day1.Lt',
	'X400.Day2.Dk',
	'X400.Day2.Lt',
	'X400.Day3.Dk',
	'X800.Day1.Dk',
	'X800.Day1.Lt',
	'X800.Day2.Dk',
	'X800.Day2.Lt',
	'X800.Day3.Dk'
)

late = c(
	'X400.Day3.Lt',
	'X400.Day4.Dk',
	'X400.Day4.Lt',
	'X400.Day5.Dk',
	'X400.Day5.Lt',
	'X800.Day3.Lt',
	'X800.Day4.Dk',
	'X800.Day4.Lt',
	'X800.Day5.Dk'
)

permutation.test =
	function(data,gids,conds1,conds2,niter=1000)
{
	genes = data[ data$geneid %in% gids, ]
	change = mean( apply( genes, 1, function(x){ mean( as.numeric(x[conds2])) - mean( as.numeric(x[conds1])) } ) )
	ngenes = nrow(genes)
	cat( sprintf( 'change for conds:\n%s\nvs.\n%s\nis:%g (%i genes)\n', paste(conds2,collapse=' '), paste(conds1,collapse=' '), change, ngenes ) )

	tests = sapply( 1:niter, function(x){
		rnd.genes = data[ sample(1:nrow(data),ngenes), ]
		mean( apply( rnd.genes, 1, function(x){ mean( as.numeric(x[conds2])) - mean( as.numeric(x[conds1])) } ) )
	})

	higher = length( which(tests >= change) ) / niter
	lower  = length( which(tests <= change) ) / niter
	cat( sprintf( '%i iters:\nfraction higher:%g\nfraction lower:%g\n', niter, higher, lower ) )
	return( list( change=change, p.higher=higher, p.lower=lower, ngenes=ngenes, niter=niter, dist=tests ) )
}

