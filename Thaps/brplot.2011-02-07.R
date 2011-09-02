### Justin Ashworth 2010
### Institute for Systems Biology

### Data processing and plotting for manual and automatic measurements of Thalassiosira bioreactor culture

dec.hours =
	function(x)
	# H:M:S to decimal hours
{
	s=strsplit(as.character(x),':')
	as.numeric(s[[1]][1]) + as.numeric(s[[1]][2])/60 + as.numeric(s[[1]][3])/60/60
}

unfactor =
	function(x){levels(x)[as.numeric(x)]}

delim.expand =
	function(x,y,delim=';')
	# expand paired data in which y has multiple delimited values for individual x values
	# works on (and filters/extends) a pair of corresponding vectors
	# simplifies bookkeeping and data entry
{
	if(is.factor(x)){x=unfactor(x)}
	if(is.factor(y)){y=unfactor(y)}
	for(i in grep(delim,y)){
		sp=unlist(strsplit(as.character(y[i]),delim))
		for(s in sp){
			#cat(x[i],s,'\n')
			x=c(x,x[i])
			y=c(y,s)
		}
	}
	xy=cbind(x,y)
	return(xy[!grepl(delim,y) & !y=="",])
}

delim.mean.sd =
	function(x,delim=';')
	# not vectorized: use an apply method
{
	sp=unlist(strsplit(as.character(x),delim))
	return(c(mean(as.numeric(sp),na.rm=T),sd(as.numeric(sp),na.rm=T)))
}

####################################################################################################

manual=read.delim('manual.data.txt')
auto=read.delim('auto.data.txt')

# reformat manual time
dt=strptime(manual$date.time,'%m/%d/%y %I:%M %p'); manual$date.time=dt; rm(dt)
# convert expt.time to decimal hours
hours=sapply(manual$expt.time,dec.hours); manual$expt.time=hours; rm(hours)

maxhours=max(manual$expt.time)-19
# truncate logged data
auto=subset(auto,LogTime<maxhours)

# shift (adjustment) between auto and manual times
timesynch=-0.25

# synchronize manual and auto data
auto$LogTime=auto$LogTime + timesynch

#xlim=c(min(manual$expt.time,auto$LogTime),max(manual$expt.time,auto$LogTime))
xlim=c(16,maxhours+timesynch)


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
	cells=c(-0,3),
	pH=c(6.4,8.6),
#	dO2=c(0,6),
	dO2=c(60,110),
	flu=c(0,0.7e4),
	PE=c(0,0.7),
	Temp=c(19.8,28),
	PO4=c(0,max( as.numeric(manual$PO4,na.rm=T ))/nutrient_plot_factor),
	SiO4=c(0,max(as.numeric(manual$SiO4,na.rm=T))/nutrient_plot_factor),
	NO3=c(0,max( as.numeric(manual$NO3,na.rm=T ))/nutrient_plot_factor),
	NH4=c(0,max( as.numeric(manual$NH4,na.rm=T ))/nutrient_plot_factor)
)

size='screen'
size='print'

if(size=='screen'){
# for 'screen size'
lwd=3
par(lwd=lwd,cex.axis=1.1,cex.main=2)
par(mar=c(5,6,4,14))
par(xaxs='i')
cex.mtext=1.5
}

if(size=='print'){
# for 'print size'
# e.g. pdf('',width=21,height=16)
lwd=6
par(lwd=lwd,cex.axis=2.2,cex.main=4)
par(mar=c(8,10,6,20))
par(xaxs='i')
cex.mtext=3
}

# expand delimited multiple values
cells.vs.time=delim.expand(manual$expt.time,manual$cells...100.nL)
# convert cells/100nL to 10^6 cells/mL
cells.vs.time[,2]=as.numeric(cells.vs.time[,2])/100
# base plot window (cell count is left y-axis)
main=expression('Response to rapid dynamic '*CO[2]*' enrichment')
plot(cells.vs.time,xlim=xlim,ylim=ylims$cells,cex.lab=cex.mtext*1.33,cex.main=cex.mtext*1.1,main=main,xlab='',ylab='',type='n',yaxt='n')
line=3
if(size=='print'){line=5}
mtext('Time (h)',side=1,line=line,cex=cex.mtext*1.33)

axis(2,pretty(c(0,ylims$cells[2])),col=colors$cells,col.axis=colors$cells,lwd=lwd)
abline(v=0,lty=2)
#text(4,2.5,'[inoculation]',cex=cex.mtext,srt=90)
line=3
if(size=='print'){line=6}
mtext(expression(10^6~'cells/mL'),side=2,line=line,col=colors$cells,cex=cex.mtext*1)

lowessf=0.1
# plot cells/mL
#points(cells.vs.time,pch=16,col=colors$cells)
# plot average trace
cells.vs.time.avg=t(sapply(manual$cells...100.nL,delim.mean.sd))
cells.vs.time.avg=cbind(manual$expt.time,cells.vs.time.avg[,1]/100)
colnames(cells.vs.time.avg)=c('time','cells')
points(cells.vs.time.avg,pch=16,col=colors$cells)
lines(lowess(na.exclude(cells.vs.time),f=lowessf),lty=2,col=colors$cells.trans)

# pH
# adjust bioreactor probe data by correspondence to manual auxiliary probe data
# bioreactor probe values are adjusted by mean difference from mCP pH values within relevant time window
mCP.vs.time.avg=t(sapply(manual$pH.mCP,delim.mean.sd))
mCP.vs.time.avg=cbind(manual$expt.time,mCP.vs.time.avg[,1])
relevant=manual$expt.time > 6 & manual$expt.time < maxhours

# simple scale shift
adj=c(-.35,1)
#adj=mean(mCP.vs.time.avg[relevant,2]-manual$pH.probe[relevant],na.rm=T)
adj=c(adj,0)
lmfit=FALSE
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

# probe line
lines(auto$LogTime,sapply(auto$pH,function(x){adj[2]*x+adj[1]}),col=colors$pH)
#points(manual$expt.time,manual$pH.probe+adj,col=colors$pH,pch=16)
padj=NA
if(size=='print'){padj=0.8}
line=1
if(size=='print'){line=2}
axis(4,seq(7.1,8.3,0.2),lwd=lwd,padj=padj,line=line)
line=1
if(size=='print'){line=4}
mtext('pH',side=4,line=line,at=8.6,cex=cex.mtext*0.8)

# plot dO2
par(new=T)
# base plot window
plot(manual$expt.time,manual$dO2,xlim=xlim,ylim=ylims$dO2,ylab='',xlab='',col=colors$dO2,axes=F,type='n',bg='transparent')

# probe reading (raw %)
lines(auto$LogTime,auto$dO2,col=colors$dO2)
# line for simple factor adjustment
#lines(auto$LogTime,auto$dO2*adj,col=colors$dO2)

line=4
if(size=='print'){line=8}
#axis(4,pretty(c(85,100)),line=line,col=colors$dO2,col.axis=colors$dO2,lwd=lwd,padj=padj)
#mtext(expression(paste(dO[2],' (%)')),4,adj=0.85,line=line+1,col=colors$dO2,cex=cex.mtext*0.8)
text(27,90,expression(dO[2]),col=colors$dO2,cex=cex.mtext*0.7)

# plot temperature
par(new=T)
plot(auto$LogTime,auto$Temp,xlim=xlim,ylim=ylims$Temp,ylab='',xlab='',col=colors$Temp,axes=F,type='l',bg='transparent')
line=1
if(size=='print'){line=2}
axis(4,c(19.5,20.5),line=line,col=colors$Temp,col.axis=colors$Temp,lwd=lwd,padj=padj)
line=1
if(size=='print'){line=4}
mtext(expression(paste({}^o,'C')),4,line=line,at=21.2,col=colors$Temp,cex=cex.mtext*0.8)

# plot fluorescence
par(new=T)
flu=t(sapply(manual$flu,delim.mean.sd))
flu.vs.time=delim.expand(manual$expt.time,manual$flu)
#plot(na.exclude(cbind(manual$expt.time,flu[,1])),xlim=xlim,ylim=ylims$flu,ylab='',xlab='',col=colors$flu.trans,axes=F,pch=16,bg='transparent')
plot(flu.vs.time,xlim=xlim,ylim=ylims$flu,ylab='',xlab='',col=colors$flu,axes=F,pch=16,bg='transparent')
#lines(lowess(na.exclude(cbind(manual$expt.time,flu[,1])),f=lowessf),lty=2,col=colors$flu)
lines(lowess(na.exclude(flu.vs.time),f=lowessf*0.8),lty=2,col=colors$flu.trans)
text(20,2000,'fluorescence',col=colors$flu,cex=cex.mtext*0.8)

# plot photoefficiency
dcmu=t(sapply(manual$dcmu,delim.mean.sd))
manual$pe=(dcmu[,1]-flu[,1])/dcmu[,1]
#manual$pe.sd=?
par(new=T)
#ylims$PE=c(0,max(manual$pe,na.rm=T))
plot(na.exclude(cbind(manual$expt.time,manual$pe)),xlim=xlim,ylim=ylims$PE,ylab='',xlab='',col=colors$PE,axes=F,pch=16,bg='transparent')
lines(lowess(na.exclude(cbind(manual$expt.time,manual$pe)),f=lowessf*1),lty=2,col=colors$PE.trans)
#axis(4,c(0.3,0.4,0.5,0.6),line=line,col=colors$PE,col.axis=colors$PE,lwd=lwd,padj=padj)
#mtext('PhotoEff.',4,adj=1.05,line=line+1,col=colors$PE,cex=cex.mtext*0.8)
text(24,0.58,'photoefficiency',col=colors$PE,cex=cex.mtext*0.8)

# mark RNA sampling
rna=subset(manual,RNA.mL.!="",select=c(expt.time,RNA.mL.))
if(!is.null(manual$RNA.arrays)){
	rna=subset(manual,RNA.arrays!="",select=c(expt.time,RNA.arrays))
}
arrows.start=ylims$PE[2]*1.05
arrows.end=ylims$PE[2]
arrows(rna$expt.time,rep(arrows.start,nrow(rna)),rna$expt.time,rep(arrows.end,nrow(rna)),length=0.03*cex.mtext,lwd=cex.mtext)
line=-3
if(size=='print'){line=-5}
mtext('RNA sampled',3,line=line,cex=cex.mtext,adj=0.6)

# note equilibration period
line=-1
if(size=='print'){line=-5}
#mtext('[equilibration]',line=line,side=2,cex=cex.mtext,adj=0.9)

nutrients=FALSE
if(nutrients){
	# nutrients
	trendlines=FALSE
	trendlines=TRUE
	trendf=0.1
	legx=xlim[2]-10

	# plot PO4
	par(new=T)
	pcex=1
	if(size=='print'){pcex=1.5}
	plot(na.exclude(cbind(manual$expt.time,manual$PO4)),xlim=xlim,ylim=ylims$PO4,ylab='',xlab='',col=colors$PO4,axes=F,pch=16,bg='transparent',cex=pcex)
	if(trendlines){lines(lowess(na.exclude(cbind(manual$expt.time,manual$PO4)),f=trendf),lty=2,col=colors$PO4)}
	text(legx,ylims$PO4[2]*0.1,expression(PO[4]),col=colors$PO4,cex=cex.mtext)
	#axis(4,pretty(c(18,22)),line=3,col=colors$PO4,col.axis=colors$PO4,lwd=lwd)
	#mtext(expression(micro~M),4,line=6,at=20,col=colors$PO4,cex=cex.mtext)

	# plot SiO4
	par(new=T)
	plot(na.exclude(cbind(manual$expt.time,manual$SiO4)),xlim=xlim,ylim=ylims$SiO4,ylab='',xlab='',col=colors$SiO4,axes=F,pch=16,bg='transparent',cex=pcex)
	if(trendlines){lines(lowess(na.exclude(cbind(manual$expt.time,manual$SiO4)),f=trendf),lty=2,col=colors$SiO4)}
	text(legx,ylims$SiO4[2]*0.15,expression(SiO[4]),col=colors$SiO4,cex=cex.mtext)
	#axis(4,pretty(c(18,22)),line=3,col=colors$SiO4,col.axis=colors$SiO4,lwd=lwd)
	#mtext(expression(micro~M),4,line=6,at=20,col=colors$SiO4,cex=cex.mtext)

	# plot NO3
	par(new=T)
	plot(na.exclude(cbind(manual$expt.time,manual$NO3)),xlim=xlim,ylim=ylims$NO3,ylab='',xlab='',col=colors$NO3,axes=F,pch=16,bg='transparent',cex=pcex)
	if(trendlines){lines(lowess(na.exclude(cbind(manual$expt.time,manual$NO3)),f=trendf),lty=2,col=colors$NO3)}
	text(legx,ylims$NO3[2]*0.2,expression(NO[3]),col=colors$NO3,cex=cex.mtext)
	#axis(4,pretty(c(18,22)),line=3,col=colors$NO3,col.axis=colors$NO3,lwd=lwd)
	#mtext(expression(micro~M),4,line=6,at=20,col=colors$NO3,cex=cex.mtext)

	# plot NH4
	par(new=T)
	plot(na.exclude(cbind(manual$expt.time,manual$NH4)),xlim=xlim,ylim=ylims$NH4,ylab='',xlab='',col=colors$NH4,axes=F,pch=16,bg='transparent',cex=pcex)
	if(trendlines){lines(lowess(na.exclude(cbind(manual$expt.time,manual$NH4)),f=trendf),lty=2,col=colors$NH4)}
	text(legx,ylims$NH4[2]*0.25,expression(NH[4]),col=colors$NH4,cex=cex.mtext)
	#axis(4,pretty(c(18,22)),line=3,col=colors$NH4,col.axis=colors$NH4,lwd=lwd)
	#mtext(expression(micro~M),4,line=6,at=20,col=colors$NH4,cex=cex.mtext)
}
