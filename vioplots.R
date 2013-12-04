require(ggplot2)
#d = read.delim('Es.2')
#names(d) = c('class','ind','score')
#x=d

vp_plot = function(x,writepdf=T,versions=10,do.t=T){
	x$class = factor(x$class, levels=rev(levels(x$class)))
	p = ggplot(x, aes(factor(class),score))
	if(do.t) tt = t.test(x$score~x$class)
	xpos = 1.5
	ypos = 0.9*range(x$score)[2]+0.1*range(x$score)[1]
#	cat(xpos, ypos,'\n')
	for(i in 1:versions){
		if(writepdf) pdf(paste('energies',i,'pdf',sep='.'),width=4,height=6,useDingbats=F)
		else dev.new()
			pp = ( p
				+ geom_violin(trim=F)
				+ geom_jitter(position=position_jitter(width=0.2),shape=16)
				)
			if(do.t) pp = pp + annotate("text", x=xpos, y=ypos, label=sprintf('p = %.3g',tt$p.value))
			print(pp)
		if(writepdf) dev.off()
	}
}
