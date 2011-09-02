source('brplot.R')

# pH bioreactor probe vs. mCP measurements
# ...

# fluorescence vs. cells/mL
flu.vs.time.avg=t(sapply(manual$fluorescence,delim.mean.sd))
flu.vs.time.avg=cbind(manual$expt.time,flu.vs.time.avg[,1])

# reduce data to precisely corresponding timepoints (probably more strict than necessary)
fluex=na.exclude(flu.vs.time.avg)
cellex=na.exclude(cells.vs.time.avg)
fluex[,1] %in% cellex[,1]
cellex[,1] %in% fluex[,1]
fluex=fluex[fluex[,1] %in% cellex[,1],]
cellex=cellex[cellex[,1] %in% fluex[,1],]
fluex[,1]==cellex[,1]
flu.vs.cells=cbind(cellex[,1],cellex[,2],fluex[,2])
fit=lm(flu.vs.cells[,3]~flu.vs.cells[,2])
# partition data into day and night
# note: 12.7 is experiment-specific!
lightsT = (flu.vs.cells[,1] %% 24) - 12.7 < 0
lightsT
png('fluorescence.vs.cells.per.mL.png',width=800,height=800)
par(lwd=2)
par(cex=2)
plot(flu.vs.cells[,2],flu.vs.cells[,3],type='n',xlab=expression(10^6~cells/mL),ylab='fluoresence (485nm->680nm)',main='Thaps fluoresence vs. cells/mL')
abline(fit)
points(flu.vs.cells[!lightsT,2],flu.vs.cells[!lightsT,3],col='blue',pch=20)
points(flu.vs.cells[lightsT,2],flu.vs.cells[lightsT,3],col='orange',pch=20)
dev.off()
