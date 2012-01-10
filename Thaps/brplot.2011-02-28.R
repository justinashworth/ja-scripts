### Justin Ashworth 2010
### Institute for Systems Biology

### Data processing and plotting for manual and automatic measurements of Thalassiosira bioreactor culture

source('brplot.util.R')

####################################################################################################

manual=read.delim('manual.data.txt')
#auto=read.delim('auto.data.txt')

# processing the time column of the large automatic data set is time-consuming--do it once, and save/reuse processed file
process_auto_time=TRUE
process_auto_time=FALSE
auto=read.delim('auto.processed.txt')

# reformat manual time
dt=strptime(manual$date.time,'%m/%d/%y %I:%M %p'); manual$date.time=dt; rm(dt)
# convert expt.time to decimal hours
hours=sapply(manual$expt.time,dec.hours); manual$expt.time=hours; rm(hours)

if(process_auto_time){
	# convert to decimal hours (note: do only if necessary)
	hours=sapply(auto$LogTime,dec.hours); auto$LogTime=hours; rm(hours)
	# filter auto data (reduce sampling interval from 1 min to ~10)
	auto=subset(auto,LogTime %% .17 < 0.02)
	write.table(auto,'auto.processed.txt',quote=F,row.names=F,sep='\t')
}

maxhours=162
# truncate logged data
auto=subset(auto,LogTime<maxhours)

# note: what are these magic numbers?
#timesynch=75.3-99.1
timesynch=-0.5

# synchronize manual and auto data
auto$LogTime=auto$LogTime + timesynch

#xlim=c(min(manual$expt.time,auto$LogTime),max(manual$expt.time,auto$LogTime))
#xlim=c(min(manual$expt.time,auto$LogTime),maxhours+timesynch)
xlim=c(-10,maxhours+timesynch-18)
#xlim=c(min(manual$expt.time,auto$LogTime),173.75+timesynch)

####################################################################################################
### PLOT
# use this list for master control of color arguments in plotting commands
colors=list(
	cells='#0000FF',
	cells.trans='#0000FF66',
	pH='black',
	pH.mCP='#55555599',
	dO2='red',
	flu='#006600',
	flu.trans='#00660066',
	PE='#00CC00',
	PE.trans='#00CC0066',
	Temp='magenta',
	night='lightgray',
	day='white',
	PO4='orange',
	SiO4='cyan',
	NO3='purple',
	NH4='darkblue'
)

nutrient_plot_factor=0.3

ylims=list(
	cells=c(-0.5,3.2),
	pH=c(6.4,9),
#	dO2=c(0,6),
	dO2=c(50,120),
	flu=c(-5e3,0.9e4),
	PE=c(-0.2,0.7),
	Temp=c(12,36),
	PO4=c(0,max( as.numeric(manual$PO4 ),na.rm=T)/nutrient_plot_factor),
	SiO4=c(0,max(as.numeric(manual$SiO4),na.rm=T)/nutrient_plot_factor),
	NO3=c(0,max( as.numeric(manual$NO3 ),na.rm=T)/nutrient_plot_factor),
	NH4=c(0,max( as.numeric(manual$NH4 ),na.rm=T)/nutrient_plot_factor)
)

mode='screen'
mode='print'

if(mode=='screen'){
# for 'screen mode'
lwd=3
par(lwd=lwd,cex.axis=1.1,cex.main=2)
par(mar=c(5,6,4,14))
par(xaxs='i')
cex.mtext=1.5
cex.points=1
}

if(mode=='print'){
# for 'print mode'
# e.g. pdf('',width=21,height=16)
lwd=6
par(lwd=lwd,cex.axis=2.2,cex.main=4)
par(mar=c(8,10,6,20))
par(xaxs='i')
cex.mtext=3
cex.points=2
}

# expand delimited multiple values
cells.vs.time=delim.expand(manual$expt.time,manual$cells...100.nL)
# convert cells/100nL to 10^6 cells/mL
cells.vs.time[,2]=as.numeric(cells.vs.time[,2])/100
# base plot window (cell count is left y-axis)
#ylims$cells=as.numeric(c(min(cells.vs.time[,2],na.rm=T),as.numeric(max(cells.vs.time[,2],na.rm=T))+topinnermargin))
#main=expression(paste(bolditalic('Thalassiosira'),' growth in elevated ',CO[2],' (800 ppm)'))
main=expression(bolditalic('T. pseudonana')*' growth in normal '*CO[2]*' (~400 ppm)')
plot(cells.vs.time,xlim=xlim,ylim=ylims$cells,cex.lab=cex.mtext*1.33,main=main,xlab='',ylab='',type='n',yaxt='n')
line=3
if(mode=='print'){line=5}
mtext('Time (h)',side=1,line=line,cex=cex.mtext*1.33)

axis(2,pretty(c(0,ylims$cells[2])),col=colors$cells,col.axis=colors$cells,lwd=lwd)
abline(v=0,lty=2)
text(4,2.8,'[inoculation]',cex=cex.mtext,srt=90)
line=3
if(mode=='print'){line=6}
mtext(expression(10^6~'cells/mL'),side=2,line=line,col=colors$cells,cex=cex.mtext*1.33)

daynight=TRUE
if(daynight){
	# plot day/night (must be plotted first)
	stop=xlim[2]
	stop=xlim[2]*1.05
	stop=maxhours
	lightsout=seq(12,stop,24)
	lightson=seq(24,stop,24)
	while(length(lightson) < length(lightsout)){
		i=length(lightson)+1
		lightson[i]=lightsout[i]
	}
	while(length(lightsout) < length(lightson)){
		i=length(lightsout)+1
		lightsout[i]=lightson[i]
	}
	#if( max(lightsout) > max(lightson) ) {
	#	lightson=c(lightson,stop)
	#} else {
	#	lightsout=c(lightsout,stop)
	#}

	# time correction (?)
	tcor=0.7
	lightson=lightson+tcor
	lightsout=lightsout+tcor
	rect(lightsout,rep(ylims$cells[1]-0.1,length(lightsout)),lightson,rep(ylims$cells[2]+0.1,length(lightson)),col=colors$night,lty=0)
}

lowessf=0.1

# plot cells/mL
#points(cells.vs.time,pch=16,col=colors$cells)
# plot average trace
cells.vs.time.avg=t(sapply(manual$cells...100.nL,delim.mean.sd))
cells.vs.time.avg=cbind(manual$expt.time,cells.vs.time.avg[,1]/100)
colnames(cells.vs.time.avg)=c('time','cells')
#cells.vs.time.avg=as.data.frame(na.exclude(cells.vs.time.avg))
points(cells.vs.time.avg,pch=16,col=colors$cells.trans,cex=cex.points)
#lines(lowess(na.exclude(cells.vs.time.avg),f=lowessf),lty=2,col=colors$cells)
lines(lowess(na.exclude(cells.vs.time),f=lowessf),lty=2,col=colors$cells)
# don't average?
#lines(lowess(na.exclude(manual$cells...100.nL),f=0.5),lty=2,col=colors$cells)
#text(xlim[2],1.2,'cells',col=colors$cells)

# pH
# adjust bioreactor probe data by correspondence to manual auxiliary probe data
# bioreactor probe values are adjusted by mean difference from mCP pH values within relevant time window
mCP.vs.time.avg=t(sapply(manual$pH.mCP,delim.mean.sd))
mCP.vs.time.avg=cbind(manual$expt.time,mCP.vs.time.avg[,1])
relevant=manual$expt.time > 6 & manual$expt.time < maxhours

adj=c(-.35,1)
# simple scale shift
adj=mean(mCP.vs.time.avg[relevant,2]-manual$pH.probe[relevant],na.rm=T)
adj=c(adj,0)
lmfit=TRUE
if(lmfit){
	# try linear fit (mx+b)
	fit=lm(mCP.vs.time.avg[relevant,2] ~ manual$pH.probe[relevant])
	adj=fit$coefficients
}
# plot pH
par(new=T)
# (bioreactor probe)
plot(auto$LogTime,auto$pH,xlim=xlim,ylim=ylims$pH,ylab='',xlab='',axes=F,type='n',bg='transparent')
# hack to restore right 'plot border' that is sometimes covered up by a gray rectangle when axes=F and par(xaxs='i')
abline(v=xlim[2])

pHbaseline=FALSE
if(pHbaseline){
	# starting media baseline
	pHbasetime=c(-30,-25)
	pHbaseline=mean(auto$pH[auto$LogTime>pHbasetime[1] & auto$LogTime<pHbasetime[2]])
	#abline(h=pHbaseline+adj,col='gray',lty=2,lwd=2)
	abline(h=pHbaseline*adj[2]+adj[1],col='gray',lty=2,lwd=2)
}

# pH (mCP assay)
points(mCP.vs.time.avg,col=colors$pH.mCP,pch=8,cex=0.1)
#mCP=delim.expand(manual$expt.time,manual$pH.mCP)
#points(mCP,col=colors$pH.mCP,pch=8,cex=0.1)

# probe line
#lines(auto$LogTime,auto$pH+adj,col=colors$pH)
lines(auto$LogTime,sapply(auto$pH,function(x){adj[2]*x+adj[1]}),col=colors$pH)
#points(manual$expt.time,manual$pH.probe+adj,col=colors$pH,pch=16)
padj=NA
if(mode=='print'){padj=0.8}
line=1
if(mode=='print'){line=2}
axis(4,pretty(c(7.4,8.4)),lwd=lwd,padj=padj,line=line)
line=2
if(mode=='print'){line=4}
mtext('pH',side=4,line=line,at=8.6,cex=cex.mtext*1)

# dO2
dO2.manual=FALSE
if(dO2.manual){
	# convert probe voltages to mg/L
	v.to.sal.slope=4.4789
	v.to.sal.inter=-1.3437
	# salinity, absolute fraction (parts-per-thousand(ppt)/1000)
	salinity=0.032
	# factor by which salinity affects dO2 reading
	sal.factor=49.64
	manual$dO2 = (manual$dO2.electrode * v.to.sal.slope) + v.to.sal.inter - (salinity * sal.factor)

	# attempt to calibrate the bioreactor probe trace to the external measurements
	# simple factor.
	adj=1
	adj=mean(manual$dO2 / manual$dO2..., na.rm=TRUE)
	# linear fit (mx+b)
	#fit=lm(manual$dO2...[relevant] ~ manual$dO2[relevant])
	fit=lm(manual$dO2[relevant] ~ manual$dO2...[relevant])
}

# plot dO2
par(new=T)
# base plot window
plot(manual$expt.time,manual$dO2,xlim=xlim,ylim=ylims$dO2,ylab='',xlab='',col=colors$dO2,axes=F,type='n',bg='transparent')

# probe reading (raw %)
lines(auto$LogTime,auto$dO2,col=colors$dO2)
# line for simple factor adjustment
#lines(auto$LogTime,auto$dO2*adj,col=colors$dO2)
# the manual probe measurements are all over the place--unsuitable for now
#lines(auto$LogTime,sapply(auto$dO2,function(x){fit$coefficients[2]*x+fit$coefficients[1]}),col=colors$dO2)

# bioreactor dO2 probe calibration (source: Elizabeth?)

#points(manual$expt.time,manual$dO2,col=colors$dO2)
line=4
if(mode=='print'){line=8}
axis(4,pretty(c(85,100)),at=pretty(c(85,100)),labels=pretty(c(85,100)),line=line,col=colors$dO2,col.axis=colors$dO2,lwd=lwd,padj=padj)
mtext(expression(paste(dO[2],' (%)')),4,adj=0.85,line=line+1,col=colors$dO2,cex=cex.mtext*0.8)
#mtext(expression(paste(dO[2],' (mg/L)')),4,adj=1,line=line,col=colors$dO2,cex=cex.mtext*0.8)

# plot temperature
plot_temp = FALSE
if ( plot_temp ) {
par(new=T)
plot(auto$LogTime,auto$Temp,xlim=xlim,ylim=ylims$Temp,ylab='',xlab='',col=colors$Temp,axes=F,type='l',bg='transparent')
line=4
if(mode=='print'){line=8}
axis(4,pretty(c(18,22)),line=8,col=colors$Temp,col.axis=colors$Temp,lwd=lwd,padj=padj)
line=7
if(mode=='print'){line=14}
mtext(expression(paste('Temp. (',{}^o,'C)')),4,line=line,at=20,col=colors$Temp,cex=cex.mtext)
}

# plot fluorescence
par(new=T)
flu=t(sapply(manual$flu,delim.mean.sd))
flu.vs.time=delim.expand(manual$expt.time,manual$flu)
#plot(na.exclude(cbind(manual$expt.time,flu[,1])),xlim=xlim,ylim=ylims$flu,ylab='',xlab='',col=colors$flu.trans,axes=F,pch=16,bg='transparent')
plot(flu.vs.time,xlim=xlim,ylim=ylims$flu,ylab='',xlab='',col=colors$flu.trans,axes=F,pch=16,cex=cex.points,bg='transparent')
#lines(lowess(na.exclude(cbind(manual$expt.time,flu[,1])),f=lowessf),lty=2,col=colors$flu)
lines(lowess(na.exclude(flu.vs.time),f=lowessf*0.8),lty=2,col=colors$flu)
text(xlim[2]-50,max(na.exclude(flu))+1000,'fluorescence',col=colors$flu,cex=cex.mtext)

# plot photoefficiency
dcmu=t(sapply(manual$dcmu,delim.mean.sd))
manual$pe=(dcmu[,1]-flu[,1])/dcmu[,1]
manual$pe[manual$pe < 0] = NA
#manual$pe.sd=?
par(new=T)
line=7
if(mode=='print'){line=14}
#ylims$PE=c(0,max(manual$pe,na.rm=T))
plot(na.exclude(cbind(manual$expt.time,manual$pe)),xlim=xlim,ylim=ylims$PE,ylab='',xlab='',col=colors$PE.trans,axes=F,pch=16,cex=cex.points,bg='transparent')
lines(lowess(na.exclude(cbind(manual$expt.time,manual$pe)),f=lowessf*1),lty=2,col=colors$PE)
axis(4,c(0.3,0.4,0.5,0.6),line=line,col=colors$PE,col.axis=colors$PE,lwd=lwd,padj=padj)
mtext('PhotoEff.',4,adj=1.05,line=line+1,col=colors$PE,cex=cex.mtext*0.8)

# mark RNA sampling
rna=subset(manual,RNA.mL.!="",select=c(expt.time,RNA.mL.))
if(!is.null(manual$RNA.arrays)){
	rna=subset(manual,RNA.arrays!="",select=c(expt.time,RNA.arrays))
}
arrows.start=ylims$PE[2]*1.05
arrows.end=ylims$PE[2]
arrows(rna$expt.time,rep(arrows.start,nrow(rna)),rna$expt.time,rep(arrows.end,nrow(rna)),length=0.03*cex.mtext,lwd=cex.mtext)
line=-2
if(mode=='print'){line=-6}
mtext('cells sampled (for RNA)',3,line=line,cex=cex.mtext)

# note equilibration period
line=-1
if(mode=='print'){line=-4}
mtext('[equilibration]',line=line,side=2,cex=cex.mtext,adj=0.95)

nutrients=FALSE
nutrients=TRUE
if(nutrients){
	# nutrients
	trendlines=FALSE
	trendlines=TRUE
	trendf=0.1
	legx=xlim[2]-15

	# plot PO4
	par(new=T)
	pcex=1
	if(mode=='print'){pcex=1.5}
	plot(na.exclude(cbind(manual$expt.time,manual$PO4)),xlim=xlim,ylim=ylims$PO4,ylab='',xlab='',col=colors$PO4,axes=F,pch=16,bg='transparent',cex=pcex)

	maxnut = max(as.numeric(manual$PO4),na.rm=TRUE)
	abline(h=maxnut,lty=2,lwd=lwd)
	text(100,maxnut*0.9,'100%',cex=cex.mtext)

	if(trendlines){lines(lowess(na.exclude(cbind(manual$expt.time,manual$PO4)),f=trendf),lty=2,col=colors$PO4)}
	nuttext = bquote( PO[4] (.(round(maxnut))) )
	text(legx,ylims$PO4[2]*0.12, nuttext,col=colors$PO4,cex=cex.mtext)
	#axis(4,pretty(c(18,22)),line=3,col=colors$PO4,col.axis=colors$PO4,lwd=lwd)
	#mtext(expression(micro~M),4,line=6,at=20,col=colors$PO4,cex=cex.mtext)

	# plot SiO4
	par(new=T)
	plot(na.exclude(cbind(manual$expt.time,manual$SiO4)),xlim=xlim,ylim=ylims$SiO4,ylab='',xlab='',col=colors$SiO4,axes=F,pch=16,bg='transparent',cex=pcex)
	if(trendlines){lines(lowess(na.exclude(cbind(manual$expt.time,manual$SiO4)),f=trendf),lty=2,col=colors$SiO4)}
	maxnut = max(as.numeric(manual$SiO4),na.rm=TRUE)
	nuttext = bquote( SiO[4] (.(round(maxnut))) )
	text(legx,ylims$SiO4[2]*0.18,nuttext,col=colors$SiO4,cex=cex.mtext)
	#axis(4,pretty(c(18,22)),line=3,col=colors$SiO4,col.axis=colors$SiO4,lwd=lwd)
	#mtext(expression(micro~M),4,line=6,at=20,col=colors$SiO4,cex=cex.mtext)

	# plot NO3
	par(new=T)
	plot(na.exclude(cbind(manual$expt.time,manual$NO3)),xlim=xlim,ylim=ylims$NO3,ylab='',xlab='',col=colors$NO3,axes=F,pch=16,bg='transparent',cex=pcex)
	if(trendlines){lines(lowess(na.exclude(cbind(manual$expt.time,manual$NO3)),f=trendf),lty=2,col=colors$NO3)}
	maxnut = max(as.numeric(manual$NO3),na.rm=TRUE)
	nuttext = bquote( NO[3] (.(round(maxnut))) )
	text(legx,ylims$NO3[2]*0.24,nuttext,col=colors$NO3,cex=cex.mtext)
	#axis(4,pretty(c(18,22)),line=3,col=colors$NO3,col.axis=colors$NO3,lwd=lwd)
	#mtext(expression(micro~M),4,line=6,at=20,col=colors$NO3,cex=cex.mtext)

	NH4=FALSE
	if(NH4){
	# plot NH4
	par(new=T)
	plot(na.exclude(cbind(manual$expt.time,manual$NH4)),xlim=xlim,ylim=ylims$NH4,ylab='',xlab='',col=colors$NH4,axes=F,pch=16,bg='transparent',cex=pcex)
	if(trendlines){lines(lowess(na.exclude(cbind(manual$expt.time,manual$NH4)),f=trendf),lty=2,col=colors$NH4)}
	maxnut = max(as.numeric(manual$NH4),na.rm=TRUE)
	nuttext = bquote( NH[4] (.(round(maxnut))) )
	text(legx,ylims$NH4[2]*0.1,nuttext,col=colors$NH4,cex=cex.mtext)
	#axis(4,pretty(c(18,22)),line=3,col=colors$NH4,col.axis=colors$NH4,lwd=lwd)
	#mtext(expression(micro~M),4,line=6,at=20,col=colors$NH4,cex=cex.mtext)
	}
}
