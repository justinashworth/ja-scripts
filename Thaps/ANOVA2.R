run2 = function(el,dl,pref){
	earlylate(el,pref)
	darklight(dl,pref)
}

run = function() {
	el = read.delim('mev.earlylate.txt')
	dl = read.delim('mev.darklight.txt')
	run2(el,dl,'all')
	eln = named.genes(el)
	dln = named.genes(dl)
	run2(eln,dln,'named')
	elk = known.genes(el)
	dlk = known.genes(dl)
	run2(elk,dlk,'known')
}

tsv = function(tbl,fname,...) {
	write.table(tbl,fname,row.names=F,quote=F,sep='\t',...)
}

#ratios=read.delim('ratios.changesfromavg')

trim.mev.table =
	function(table)
{
}

known.genes =
	function(ratios)
{
	ratios = ratios[ !grepl('hypothetical protein',ratios$protein) & ! ratios$protein=='' & !is.na(ratios$protein), ]
	return(ratios)
}

named.genes =
	function(ratios)
{
	ratios = ratios[ ratios$name != '' & !is.na(ratios$name), ]
	return(ratios)
}

earlylate =
	function(ratios,filepref=NULL)
{
	ratios = ratios[ order(ratios$F.ratio..early_late.,decreasing=TRUE), ]

	early = c('X400.Day1.Dk','X400.Day1.Lt','X400.Day2.Dk','X400.Day2.Lt','X400.Day3.Dk','X800.Day1.Dk','X800.Day1.Lt','X800.Day2.Dk','X800.Day2.Lt','X800.Day3.Dk')
	late = c('X400.Day3.Lt','X400.Day4.Dk','X400.Day4.Lt','X400.Day5.Dk','X400.Day5.Lt','X800.Day3.Lt','X800.Day4.Dk','X800.Day4.Lt','X800.Day5.Dk')
	early = names(ratios) %in% early
	late = names(ratios) %in% late

	ratios$change = apply( ratios, 1, function(x) { mean(as.numeric(x[late])) - mean(as.numeric(x[early])) } )
	# reorder so ratios remain in last column positions
	ratios = ratios[ , c(which(!grepl('Dk|Lt',names(ratios))),grep('Dk|Lt',names(ratios))) ]

	upearly = ratios[ ratios$change < 0, ]
	upearly = upearly[ order(upearly$change), ]

	uplate = ratios[ ratios$change > 0, ]
	uplate = uplate[ order(uplate$change,decreasing=TRUE), ]

	if ( !is.null(filepref) ) {
		tsv( upearly, sprintf('%s.upearly',filepref) )
		tsv( upearly[1:20,], sprintf('%s.upearly.20',filepref) )
		tsv( uplate, sprintf('%s.uplate',filepref) )
		tsv( uplate[1:20,], sprintf('%s.uplate.20',filepref) )
	}

	return( list( 'upearly'=upearly, 'uplate'=uplate ) )
}

darklight =
	function(ratios,filepref=NULL)
{
	ratios = ratios[ order(ratios$F.ratio..dark_light.,decreasing=TRUE), ]
	dark = grepl('Dk',names(ratios))
	light = grepl('Lt',names(ratios))
	ratios$change = apply( ratios, 1, function(x) { mean(as.numeric(x[light])) - mean(as.numeric(x[dark])) } )
	# reorder so ratios remain in last column positions
	ratios = ratios[ , c(which(!grepl('Dk|Lt',names(ratios))),grep('Dk|Lt',names(ratios))) ]

	updark = ratios[ ratios$change < 0, ]
	updark = updark[ order(updark$change), ]

	uplight = ratios[ ratios$change > 0, ]
	uplight = uplight[ order(uplight$change,decreasing=TRUE), ]

	if ( !is.null(filepref) ) {
		tsv( updark, sprintf('%s.updark',filepref) )
		tsv( updark[1:20,], sprintf('%s.updark.20',filepref) )
		tsv( uplight, sprintf('%s.uplight',filepref) )
		tsv( uplight[1:20,], sprintf('%s.uplight.20',filepref) )
	}

	return( list( 'updark'=updark, 'uplight'=uplight ) )
}

interacting =
	function(ratios,filepref=NULL)
{
	ratios = ratios[ order(ratios$F.ratio..interaction.,decreasing=TRUE), ]
	if ( !is.null(filepref) ) {
		tsv(ratios,paste(filepref,'interacting',sep='.'))
	}
	return(ratios)
}
