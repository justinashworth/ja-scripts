source('brplot.R')
l=lowess(na.exclude(cells.vs.time.avg),f=0.1)
g=growthrate(l$x,l$y)
endpoint=length(l$x)-1
plot(l$x[1:endpoint],g)
g=g*24

png('growth.rate.png',width=800,height=600)
par(cex=2)
par(lwd=2)
plot(l$x[1:endpoint],g,ylim=c(-1,1.8),xlab='time (hrs)',ylab=expression(10^6~cells/mL/day),main="growth rate",type='n')
rect(lightsout,rep(-1.1,length(lightsout)),lightson,rep(2,length(lightson)),col=colors$night,lty=0)
points(l$x[1:endpoint],g)
lines(lowess(l$x[1:endpoint],g,f=0.1))
abline(h=0,lty=2)
dev.off()

png('doubling.time.png',width=800,height=600)
par(cex=2)
par(lwd=2)
d=log(2)/log(1+g)
dt=cbind(l$x[1:endpoint],d)
plot(dt,ylim=c(0,5),xlab='time (hrs)',ylab='doubling time (days)',type='n')
rect(lightsout,rep(-0.1,length(lightsout)),lightson,rep(5,length(lightson)),col=colors$night,lty=0)
points(dt)
final=endpoint
lines(lowess(na.exclude(dt[1:final]),f=0.4))
dev.off()
