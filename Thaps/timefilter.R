d=read.delim('alee_2010_10_11.txt')
d$hours=sapply(d$LogTime,
	function(x){
		s=strsplit(as.character(x),':')
		as.numeric(s[[1]][1]) + as.numeric(s[[1]][2])/60 + as.numeric(s[[1]][3])/60/60
	}
)
dd=subset(d,hours %% 0.2 < 0.02)
nrow(dd)
plot(dd$hours,dd$pH,pch=20)
dd$hours=round(dd$hours,2)
write.table(dd,'alee_2010_10_11.lowdens.txt',quote=FALSE,row.names=F,sep='\t')
