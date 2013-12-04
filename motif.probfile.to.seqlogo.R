# input example: position base probability matrix
##comments
#pos	A	C	G	T
#7	0	0.84	0.04	0.08
#6	0	0.94	0	0.06

# output: sequence logo

#files = readLines('prob.motifs')
#files = readLines('motif.probs.files')
files = readLines('probfiles')

require(seqLogo)

for(f in files){
	p=read.delim(f,comment.char='#')
	horiz=t(p[,c(2:ncol(p))])
	norm=t(t(horiz)*1/colSums(horiz))
	pdf(paste(f,'.pdf',sep=''),width=8,height=5)
#	seqLogo(norm)
	# turn off scaling to visually de-emphasize oversensitivity to absolute scale in the data
	#seqLogo(norm,ic.scale=F)
	seqLogo(norm,ic.scale=F)
	dev.off()
}
