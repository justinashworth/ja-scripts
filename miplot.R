# quick script to plot residue covariation from Fodor's mutual information program
fname='lrps.mi'
d=read.delim(fname)
# based on the alignment used, adjust residue numbering to match FL11 structure
d$i=d$i+25
d$j=d$j+25
# normalize scores to [0,1]
d$score=(d$score-min(d$score))/(max(d$score)-min(d$score))
# filter for top 5% covariant edges, create color gradient
d$col=d$score
e=ecdf(d$col)
plot(e)
e(0.715)
# filter
d$col[d$col<0.715]=NA
# renormalize top scores to [0,1]
d$col=(d$col-min(d$col,na.rm=T))/(max(d$col,na.rm=T)-min(d$col,na.rm=T))
# remove NA's
d$col[is.na(d$col)]=0
# replace with color gradient
d$col=rgb(1,0,0,d$col)
# plot
m=min(d$i,d$j)
x=max(d$i,d$j)
xlim=c(m,x)
ylim=xlim
m=m-10
x=x+10
plot(d$i,d$j,type='n',xlab='',ylab='',main=paste('residue pair mutual information',fname,sep='\n'),xlim=xlim,ylim=ylim)
lines(c(m,x),c(m,x))
abline(v=m:x,lty=3,col='gray')
abline(h=m:x,lty=3,col='gray')
sq=seq(m,x,by=5)
abline(h=sq,col='gray')
abline(v=sq,col='gray')
points(d$i,d$j,col=d$col,pch=15)
points(d$j,d$i,col=d$col,pch=15)
