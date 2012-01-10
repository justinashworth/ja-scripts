# Justin Ashworth
# Institute for Systems Biology

# Uses limma utilities for array normalization, maps/averages probes, calculates expression ratio changes from Day1, writes ratios file

# 1. background correction
# 2. intra-array normalization : normalizes green and red intensities
# 3. inter-array normalization : computes ratios for each sample (red) to uniform internal standard (green), and normalizes across arrays
# 4. add probe->gene mappings
# 5. average multiple values per gene
# 6. express ratios as change from Day1 within each experiment

library(limma)
# Targets.txt contains the filenames and sample ids
targets=readTargets('Targets.txt')
RG=read.maimages(targets,source='agilent')
png('array.densities.png',width=1000,height=1000)
op=par(mfrow=c(2,2))
plotDensities(RG)
RGb=backgroundCorrect(RG,method='normexp')
plotDensities(RGb)
RGn=normalizeWithinArrays(RGb,method='loess')
plotDensities(RGn)
RGn=normalizeBetweenArrays(RGn,method='Aquantile')
plotDensities(RGn)
dev.off()
# now we have intra-array, inter-array normalized probe ratios of Cy5 (red) samples to a common, uniform, independent Cy3 (green) reference (which in this case we consider an internal standard to account for array variability, rather than an experimental comparison)

# now add external probe->gene mapping (verified using BLAT) to RGn$genes table
d=read.delim('../probe.mapping')
d=d[,c(1,6)]
# this probe mapping has multiple probe->gene mappings (some probes BLAT to multiple genes)
# don't filter out probes in genes table--could be incompatible with limma subsequent methods
# instead, just map multi-mapped probes to their own probeids
tt=table(d$probe)
length(tt[tt>1])
mm=names(tt[tt>1])
d$geneid = as.character( d$geneid )
d$geneid[ d$probe %in% mm ] = as.character( d$probe[ d$probe %in% mm ] )
d = unique(d)

m=merge(RGn$genes,d,by.x='ProbeName',by.y='probe',all.x=T)
m$GeneName=as.character(m$GeneName)
m$GeneName[!is.na(m$geneid)] = as.character(m$geneid[!is.na(m$geneid)])
head(m[!is.na(m$geneid),])
m$geneid=NULL
m=m[,c(2:7,1,8:10)]
RGn$genes=m

# now average probe ratios down to gene ratio matrix
avg=avereps(RGn,ID=RGn$genes$GeneName)
ratios=avg$M

# now, subtract Day1 ratios from Day[234] ratios
# this cancels out comparison to independent Cy3 reference, and yields internally standardized ratios to Day1
# this /should/ be better than ignoring the Cy3 channel and comparing just the Cy5 (red) sample intensities directly

colnames(ratios)=targets$Sample
exp1 = c("X0329.A.Day1","X0329.A.Day2","X0329.A.Day3","X0329.A.Day4")
exp2 = c("X0329.B.Day1","X0329.B.Day2","X0329.B.Day3","X0329.B.Day4")
exp3 = c("X0922.A.Day1","X0922.A.Day2","X0922.A.Day3","X0922.A.Day4")
exp4 = c("X0922.B.Day1","X0922.B.Day2","X0922.B.Day3","X0922.B.Day4")
ratios[,exp1] = apply(ratios[,exp1],2,function(x){x-ratios[,exp1[1]]})
ratios[,exp2] = apply(ratios[,exp2],2,function(x){x-ratios[,exp2[1]]})
ratios[,exp3] = apply(ratios[,exp3],2,function(x){x-ratios[,exp3[1]]})
ratios[,exp4] = apply(ratios[,exp4],2,function(x){x-ratios[,exp4[1]]})
ratios = round(ratios,5)

# re-order columns (arrays) into chronological replicates
ratios = ratios[ , unlist( sapply( 1:4,function(x){ list(exp1[x],exp2[x],exp3[x],exp4[x]) } ) ) ]

write.table(ratios,'limma.norm.changes.from.Day1s.tsv',sep='\t')
