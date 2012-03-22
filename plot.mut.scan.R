d=read.delim('results.substitutions')
dd =d[ d$is_native==1, c('prot_chain','prot_position','aa') ]
names(dd)[ncol(dd)]='nat'
m = merge(d,dd,all.x=T,sort=F)
#m$natpos = paste( m$prot_chain, m$prot_position, m$nat)
m$natpos = paste( m$prot_position, m$nat)

chains = levels(factor(d$prot_chain))

# plot of predictions by position
for(chain in levels(factor(d$prot_chain))){
	data = m[ m$prot_chain == chain ,]

	pdf(paste('P53.3KMD.',chain,'.SNPSCAN.pdf',sep=''),width=10,height=30)
	par(cex=2,lwd=2)
	par(mfrow=c(1,2),oma=c(3,4,0,0))

	nlvls=nlevels(factor(data$natpos))
	ord = nlvls:1
	xlim=c(6,nlvls-6)
	#factors = factor(data$natpos,ordered=F)
	factors = factor(data$natpos,levels=unique(m$natpos))
	boxplot(data$bound_delta ~ factors, xlim=xlim, ylim=c(-5,40), xaxs='i', yaxs='i', las=2,pch=20,xlab='change in folding energy', horizontal=T, at=ord, names=rep('',length(ord)) )
	boxplot(data$binding_delta ~ factors, xlim=xlim, ylim=c(-5,10), xaxs='i', yaxs='i',las=2,pch=20, xlab='change in DNA binding energy', horizontal=T,at=ord )

	dev.off()
}

# inter-chain variation in predictions

# effects of missense mutations by amino acid type
muts = d[ d$is_native==0, ]

# consequences by missense type
pdf('P53.3KMD.missense.types.pdf',width=10,height=10)
par(mfrow=c(1,2),oma=c(3,4,0,0))
ord = nlevels(muts$aa):1
main='Effect of mutations TO a.a. types'
boxplot(muts$bound_delta ~ muts$aa, ylim=c(-10,30), las=2,pch=20,xlab='change in folding energy', horizontal=T, at=ord, main=main )
boxplot(muts$binding_delta ~ muts$aa, las=2,pch=20, xlab='change in DNA binding energy', horizontal=T,at=ord )
dev.off()

# consequences by native (original) type
pdf('P53.3KMD.original.types.pdf',width=10,height=10)
par(mfrow=c(1,2),oma=c(3,4,0,0))
ord = nlevels(m$nat):1
main='Effect of mutations FROM a.a. types'
boxplot(m$bound_delta ~ m$nat, las=2,pch=20,xlab='change in folding energy', horizontal=T, at=ord, main=main )
boxplot(m$binding_delta ~ m$nat, las=2,pch=20, xlab='change in DNA binding energy', horizontal=T,at=ord, names = paste('(',table(m$nat[m$is_native==1]),') ',levels(m$nat),sep='')  )
dev.off()
