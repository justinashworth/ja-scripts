library(gplots)
mt = read.table('all.mat')
mt$id = rownames(mt)
pth = read.delim('kegg.pathways')
m = merge(mt,pth)
mt = as.matrix(m[,2:5])
rownames(mt) = m$name
sig = apply(mt,1,function(x){any(x<0.05)})
sig = mt[sig,]
ord = apply(sig,1,min)
sig = sig[ order(ord), ]
lsig = log(sig,10)
minv = min(lsig)
cut = log(0.05,10)
cols = c( colorpanel(as.integer(-10*(minv-cut)),'red','yellow'), rep('gray',as.integer(-10*cut) ) )

pdf('kegg.pdf',width=17,height=22)
op=par(oma=c(4,2,4,50),cex=2,cex.lab=2,cex.axis=2,cex.main=2,cex.sub=2)
hm = heatmap.2(lsig,trace='none',scale='none',Colv=NA,Rowv=NA,col=cols,symkey=F,colsep=1:4,rowsep=1:nrow(lsig),lhei=c(1,5),lwid=c(4,3),density.info='none',cexRow=3 )
dev.off()
