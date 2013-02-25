require('seqLogo')

motif_pdfs = function(listfile)
{
	mfiles = readLines(listfile)
	cat(mfiles,'\n')
	mots = lapply(mfiles, function(x){ m = read.table(x,header=T,sep='\t'); m = m[,2:ncol(m)]; return(m) })
	names(mots) = mfiles

	for(nm in names(mots)){
		cat(nm,'\n')
		pdf(paste(nm,'pdf',sep='.'),width=8,height=4)
		seqLogo(t(mots[[nm]]))
		dev.off()
	}

	invisible(mots)
}
