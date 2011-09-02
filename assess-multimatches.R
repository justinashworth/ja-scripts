# blat data in R-readable form
d=read.delim("all.0.filtered.dat")
f=factor(d$Qname)
t=table(f)
# histogram of probe occurences in data
#hist(t,ylim=c(0,1000))
# number of probes matching multiple locations
print(length(t[which(t>1)]))

# positional occurences (matches per probe, ditributed by row order)
#png("probe-matches.png",width=1000,height=400)
#plot(t)
#dev.off()

# subsets d to remove multi-matched probes
filtered=d[which(t<2),]
chr1=subset(filtered,Tname=="Chr1")
# compute,plot resulting unprobed gaps
source("jumps.R")
j=jumps(chr1$Tstart[which(chr1$strand=="+")])
png("ok-to-remove-repeat-probes.png",width=640,height=640)
plot(ecdf(j),xlim=c(0,1000),ylim=c(0.95,1),xlab='length of unprobed region(bp)',ylab='ecdf'
,main='60-mer probes, filtered to remove probes w/ multiple matches',cex=2,cex.axis=2)
dev.off()
