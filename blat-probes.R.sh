#!/bin/sh
# $1 is
# $2 is a fasta file of appropriately named probe sequences
# $3 is
echo " \
source('~/scripts/isb-scripts/blat-probes.R')
d=read.delim('$1')
png('$3.plots.%03d.png',width=800,height=800)

# note redundant probes
#probeseqs=read.delim('$2')
require(Biostrings)
probeseqs=read.DNAStringSet('$2')
f=factor(probeseqs\$seq)
t=table(f)
redundant=t[which(t>1)]
print(paste('rendundant probes:',length(redundant)))

dnr=analMultiMatch(d,1)
write.table(dnr,'$3.nonred.dat',row.names=FALSE,quote=FALSE,sep='\t')

# write list of filtered probes
nr = probeseqs[ names(probeseqs) %in% dnr\$Qname ]
#nr = subset(probeseqs,names %in% dnr\$Qname)
write.table(nr,'$3.nr',row.names=FALSE,quote=FALSE,sep='\t')

# compute and plot gaps resulting from filtering of redundant probes
#gaps=probeGaps(d)
#nrgaps=probeGaps(dnr)
#plot(ecdf(gaps),xlim=c(0,2000),ylim=c(0.5,1),xlab='length of unprobed region(bp)',ylab='ecdf',main='distance between probe starts',cex=2,cex.axis=2)
#lines(ecdf(nrgaps),cex=2,col='green')

dev.off()
" > $3.R

R --vanilla < $3.R >> $3.log
