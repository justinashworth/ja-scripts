# for processing probe values for T. pseudonana arrays
# merges mean 3' probe values with pre-computed probe->gene map(s), arranges desriptive table including NCBI info

# processed data from Bruz's pipeline
#values=read.delim('matrix_output.tsv')
#fname='merged.intensities'
#fname='ratios.manual'
#fname='probe.ratios.tsv'
#fname='probe.ratios'
fname='probe.ratios'
cat('reading probe values...\n')
probe.values=read.delim(fname,as.is=T)

#values$probe=values$GENE
#values$GENE=NULL
#values$DESCRIPT=NULL

# strip out random and Agilent control probes
probe.values=probe.values[!grepl('^rand',probe.values$probe),]
probe.values=probe.values[!grepl('^C',probe.values$probe),]

#source('~/.Rgaggle')
#ds=getDatasetDescription()
#mat=subset(probe.values,select=c('sequence','strand','start','end',"Day1.AM_v_Day3.PM.REF.value","Day1.PM_v_Day3.PM.REF.value","Day2.AM_v_Day3.PM.REF.value","Day2.PM_v_Day3.PM.REF.value","Day3.AM_v_Day3.PM.REF.value","Day3.PM_v_Day3.PM.REF.value","Day4.AM_v_Day3.PM.REF.value","Day4.PM_v_Day3.PM.REF.value","Day5.AM_v_Day3.PM.REF.value"))
#addTrack(ds,mat,type='quantitative.segment.matrix',name='Thaps.arrays')

#chr=read.delim('probes.exons.3prime')
#org=read.delim('probes.organelles.3prime')
#probe.mapping=rbind(chr,org)
#probe.mapping=read.delim('probes.mapped')
cat('reading probe mapping...\n')
probe.mapping=read.delim('probe.mapping')
probe.mapping=probe.mapping[,c('probe','geneid')]

cat('merging probe values and mapping...\n')
mapped=merge(probe.values,probe.mapping,by=c('probe'))
ratiocolregex='ref|Day|ratio|log'
#ratiocolregex='CO2'
# slice down to geneid and values
mapped=mapped[ , c( which(names(mapped)=='geneid'), grep(ratiocolregex,names(mapped)) ) ]

#write.table(mapped,'3prime.genes.values.tsv',quote=F,row.names=F,sep='\t')

# compute average value per geneid over all probes that map to each geneid
cat('averaging probe values ~ geneid...\n')
genes = factor(mapped$geneid)
means = apply( mapped[,2:ncol(mapped)], 2, function(x) { tapply(x,factor(genes),mean) } )
means = round(means,5)
means = as.data.frame(means)
means$geneid = levels(genes)
#means = as.data.frame( cbind( levels(genes), apply( mapped[,2:ncol(mapped)], 2, function(x) { tapply(x,factor(genes),mean) } ) ) )

# unique/tapply don't seem to condense all geneids?
#means=unique(means)

#write.table(means,'3prime.genes.mean.values.tsv',quote=F,row.names=F,sep='\t')
#mat=as.matrix(means[2:ncol(means)])
#rownames(mat)=means$geneid
#source('~/.Rgaggle')
#broadcast(mat[,c(1:6)],"Thaps")
#broadcast(mat,"Thaps")

info = read.delim('tps.all.models.ja.tsv',as.is=T)
#info$id=as.character(info$gid)
#info$id[!info$genome %in% c('chloroplast','mitochondria')]=as.character(info$tag[!info$genome %in% c('chloroplast','mitochondria')])
#info$id=gsub('THAPSDRAFT_','',info$id)
#m=merge(means,info,by.x='geneid',by.y='gid')
m=merge(means,info,by.x='geneid',by.y='tag')

#m$desc=as.character(m$protein)
#m$name = as.character(m$name)
##m$name[grep('rnaseq|bd.',m$geneid)] = as.character(m$geneid[grep('rnaseq|bd.',m$geneid)])
#m$name[grep('rnaseq|bd.',m$tag)] = as.character(m$tag[grep('rnaseq|bd.',m$tag)])
#m$name[is.na(m$name)] = ''
#m$desc[m$name!='']=paste('(',m$name[m$name!=''],') ',m$desc[m$name!=''],sep='')
#
#m$desc[m$genome=='chloroplast']  = paste(m$desc[m$genome=='chloroplast'],' [chloroplast]',sep='')
m$desc[m$genome=='mitochondria'] = paste(m$desc[m$genome=='mitochondria'],' [mitochondria]',sep='')

m=m[,c(which(!grepl(ratiocolregex,names(m))),grep(ratiocolregex,names(m)))]
write.table(m,paste(fname,'mean.tsv',sep='.'),quote=F,row.names=F,sep='\t')
mm=m[ !is.na(m$protein) & m$protein != 'hypothetical protein', ]
write.table(mm,paste(fname,'mean.known.tsv',sep='.'),quote=F,row.names=F,sep='\t')
