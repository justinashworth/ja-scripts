# for processing probe ratios for T. pseudonana arrays
# merges mean 3' probe ratios with pre-computed probe->gene map(s), arranges desriptive table including NCBI info

# processed data from Bruz's pipeline
#ratios=read.delim('matrix_output.tsv')
probe.ratios=read.delim('probe.ratios.tsv')

#ratios$probe=ratios$GENE
#ratios$GENE=NULL
#ratios$DESCRIPT=NULL

# strip out random and Agilent control probes
probe.ratios=probe.ratios[!grepl('^rand',probe.ratios$probe),]
probe.ratios=probe.ratios[!grepl('^C',probe.ratios$probe),]

#source('~/.Rgaggle')
#ds=getDatasetDescription()
#mat=subset(probe.ratios,select=c('sequence','strand','start','end',"Day1.AM_v_Day3.PM.REF.ratio","Day1.PM_v_Day3.PM.REF.ratio","Day2.AM_v_Day3.PM.REF.ratio","Day2.PM_v_Day3.PM.REF.ratio","Day3.AM_v_Day3.PM.REF.ratio","Day3.PM_v_Day3.PM.REF.ratio","Day4.AM_v_Day3.PM.REF.ratio","Day4.PM_v_Day3.PM.REF.ratio","Day5.AM_v_Day3.PM.REF.ratio"))
#addTrack(ds,mat,type='quantitative.segment.matrix',name='Thaps.arrays')

chr=read.delim('probes.exons.3prime')
org=read.delim('probes.organelles.3prime')
probe.mapping=rbind(chr,org)
probe.mapping=probe.mapping[,c(1,6)]

mapped=merge(probe.ratios,probe.mapping,by=c('probe'))
# slice down to geneid and ratios
mapped=mapped[,c(which(names(mapped)=='geneid'),grep('Day',names(mapped)))]

#write.table(mapped,'3prime.genes.ratios.tsv',quote=F,row.names=F,sep='\t')

# compute average ratio per geneid over all probes that map to each geneid
means=data.frame(geneid=levels(factor(mapped$geneid)))
for(colname in grep('Day',names(mapped),value=TRUE)){
	means[[colname]]=round(tapply(mapped[[colname]],mapped$geneid,mean),5)
}

#write.table(means,'3prime.genes.mean.ratios.tsv',quote=F,row.names=F,sep='\t')
#mat=as.matrix(means[2:ncol(means)])
#rownames(mat)=means$geneid
#source('~/.Rgaggle')
#broadcast(mat[,c(1:6)],"Thaps")
#broadcast(mat,"Thaps")

info=read.delim('Thaps.genes.NCBI.tsv')
#info$id=as.character(info$gid)
#info$id[!info$genome %in% c('chloroplast','mitochondria')]=as.character(info$tag[!info$genome %in% c('chloroplast','mitochondria')])
#info$id=gsub('THAPSDRAFT_','',info$id)
m=merge(means,info,by.x='geneid',by.y='gid')

m$desc=''
m$desc[m$name!='']=paste('(',m$name[m$name!=''],') ',m$desc[m$name!=''],sep='')
m$desc[m$name!='']=paste('(',m$name[m$name!=''],') ',m$desc[m$name!=''],sep='')
m$desc[m$genome=='chloroplast']=paste(m$desc[m$genome=='chloroplast'],' [chloroplast]',sep='')
m$desc[m$genome=='mitochondria']=paste(m$desc[m$genome=='mitochondria'],' [mitochondria]',sep='')
m=m[,c(which(!grepl('Day',names(m))),grep('Day',names(m)))]
write.table(m,'ratios.3prime.mean.tsv',quote=F,row.names=F,sep='\t')
mm=m[m$protein!='hypothetical protein',]
write.table(mm,'ratios.3prime.mean.known.tsv',quote=F,row.names=F,sep='\t')
