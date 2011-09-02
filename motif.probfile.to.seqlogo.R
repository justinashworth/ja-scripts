# input example: position base probability matrix
##comments
#pos	A	C	G	T
#7	0	0.84	0.04	0.08
#6	0	0.94	0	0.06

# output: sequence logo

files = readLines('prob.motifs')

require(seqLogo)
png=1
pdf=0
if(pdf){pdf(paste(f,'.pdf',sep=''),onefile=T,width=800,height=320)}

for(f in files){
	p=read.delim(f,comment.char='#')
	horiz=t(p[,c(2:ncol(p))])
	norm=t(t(horiz)*1/colSums(horiz))
	if(png){png(paste(f,'.png',sep=''),width=800,height=320)}
#	seqLogo(norm)
	# turn off scaling to visually de-emphasize oversensitivity to absolute scale in the data
	seqLogo(norm,ic.scale=F)
	if(png){dev.off()}
}

if(pdf){dev.off()}
