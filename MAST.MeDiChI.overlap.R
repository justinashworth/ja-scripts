# compare MAST hits to Dave's MeDiChI fits on FY's new Lrp MeDiChI (ChIP) data

# for pairwise distance funcs:
source('peak.overlap.R')

# bootstrap p-value for MeDiChI fits
MeDiChI.pval.cutoff=0.05
#MeDiChI.pval.cutoff=0.02
# MAST p-value cutoff (default MAST cutoff is 1e-4, which includes some fairly weak matches)
#MAST.pval.cutoff=1e-4
# cutoff of 5e-5 usually cuts out around half of genome-wide hits
#MAST.pval.cutoff=5e-5
# don't apply cutoff here, assume MAST cutoff was sufficient
MAST.pval.cutoff=1

# keys: MAST key, values MAST and MeDiChI names
synon=list(
	'AsnC'=c('AsnC','asnC'),
	'Trh2'=c('Trh2','trh2'),
	'Trh3'=c('Trh3','trh3'),
	'Trh4'=c('Trh4','trh4'),
	'Trh6'=c('Trh6','trh6'),
	'Trh7'=c('Trh7','trh7'),
	'VNG1179C'=c('VNG1179C','p1179'),
	'VNG1237C'=c('VNG1237C','p1237')
)

conv.tf.name =
	function(x)
{
	for(key in names(synon)){
		if(x %in% synon[[key]]){
			x=key
			break
		}
	}
	return(x)
}

# keys: MAST key, values MAST and MeDiChI names
seqsyn=list(
	'Chr'=c('Chr','HALCHR'),
	'pNRC100'=c('pNRC100'),
	'pNRC200'=c('pNRC200'),
	'NC'=c('NC')
)

seqlens=list(
	'Chr'=2014239,
	'pNRC100'=191346,
	'pNRC200'=365425,
	'NC'=0
)

# R data structure for new fits file (2011-02-21_fitslist.combined.RData)
# fitslist.combined$p1179$p1179.1.1_2.1$fits.fin$NC[[1]]$coeffs.w.p.values
# position intensity p.value <
#[1,]  1442.00 0.9761781 0.5674865
#[2,]  1582.00 0.1124181 0.9784043
#[3,]  1662.00 1.1816578 0.4175165
if(!exists('MeDiChI.fits')){
	fits='2011-02-21_fitslist.combined.RData'
	cat('reading',fits,'\n')
	load(fits)
	MeDiChI.fits=fitslist.combined
}

# MAST hits flat file format
#sequence_name	motif	dir	hit_start	hit_end	score	hit_p.value
#Chr	B.sub.LrpC	fwd	14418	14432	4.25	1.09e-05
#Chr	B.sub.LrpC	rvs	17654	17668	-0.45	4.41e-05
#Chr	B.sub.LrpC	rvs	20967	20981	1.85	2.2e-05
if(!exists('masthits')){
	hits='hitsflat.tsv'
	cat('reading',hits,'\n')
	masthits=read.delim(hits)
	masthits=masthits[masthits$hit_p.value < MAST.pval.cutoff,]
}
#	masthits=masthits[grepl('Halo',masthits$motif),]

# regions in which to count overlaps
if(!exists('regions')){
	upstream=500
	downstream=100
	coords=read.delim('halo.gene.coords.tsv')
	coords$Start = as.numeric(coords$Start)
	coords$cis.start = 0
	coords$cis.end = 0
	fwd = coords$Orientation == 'For'
	coords$cis.start[fwd] = coords$Start[fwd] - upstream
	coords$cis.end[fwd] = coords$Start[fwd] + downstream
	# in halo gene coords, 'Start' (generally) refers to actual start (i.e. start > end)
	rvs = coords$Orientation == 'Rev'
	coords$cis.start[rvs] = coords$Start[rvs] + upstream
	coords$cis.end[rvs] = coords$Start[rvs] - downstream
	regions=coords[,names(coords) %in% c('canonical_Name','where','cis.start','cis.end')]
	names(regions)=c('name','seq','start','end')
}

#annotate.regions =
#	# SLOOOOW (because R sucks at for loops)
#	function(positions,regions)
#{
#	result=
#	sapply(positions,
#		function(x){
#			for(i in 1:nrow(regions)){
#				if(x>regions$start[i] & x<regions$end[i]){
#					return(as.character(regions$name[i]))
#				}
#			}
#			return(NA)
#		}
#	)
#	return(unlist(result))
#}

#if(FALSE){
##if(TRUE){
## assign regions to hits
## [got too frustrated with sapply's returning lists with this one]
#cat('assigning regions to MAST hits...\n')
#masthits$region=NA
#for(seq in names(seqlens)){
#	cat(seq,'\n')
#	seq.regions = regions[regions$seq==seq,]
#	seq.hits.i = which(masthits$sequence_name==seq)
#	seq.hitpos = (masthits$hit_start[seq.hits.i]+masthits$hit_end[seq.hits.i])/2
#	masthits$region[seq.hits.i] = annotate.regions(seq.hitpos,seq.regions)
#}}
#
#if(FALSE){
##if(TRUE){
## assign regions to peaks
#cat('assigning regions to ChIP peaks...\n')
## pre-annotate ChIP peaks in certain regions (e.g. promoters)
#for(tf in names(MeDiChI.fits)){
#	cat('tf',tf,'\n')
#	for(exp in names(MeDiChI.fits[[tf]])){
#		cat('\texp',exp,'\n')
#		seqs=names(MeDiChI.fits[[tf]][[exp]]$fits.fin)
#		for(seq in seqs){
#			if(seq=='NC')next
#			# translate to common sequence name
#			commonseq=NULL
#			for(key in names(seqsyn)){
#				if(seq %in% seqsyn[[key]]){
#					commonseq=key
#					break
#				}
#			}
#			cat('\t\tseq',commonseq,'\n')
#			seq.regions = regions[regions$seq==commonseq,]
#			peaks = MeDiChI.fits[[tf]][[exp]]$fits.fin[[seq]][[1]]$coeffs.w.p.values
#			# also apply p-value cutoff here, to avoid wasting time
#			peaks = as.data.frame(peaks[peaks[,3]<MeDiChI.pval.cutoff,])
##			print(peaks)
#			# stupid R can't do peaks[,1] on a one-row slice...
#			if(class(peaks)=='numeric') positions= round(peaks[1])
#			else positions=round(peaks[,1])
#			peak.regions=annotate.regions(positions,seq.regions)
#			MeDiChI.fits[[tf]][[exp]]$fits.fin[[seq]][[1]]$coeffs.w.p.values = cbind(peaks,peak.regions)
#		}
#	}
#}}

# count number of sites within x bp of each other
count_overlaps.distance =
	function(a,b,dis=50)
{
	distances=pairwise.site.distance(a,b)
	return(length(which(distances$dis<dis)))
}

empirical_distrib.distance =
	function(a,b,l,dis=50)
{
	if(l<1) return(0)
	distances=pairwise.site.distance( runif(a,1,l), runif(b,1,l) )
	return(length(which(distances$dis<dis)))
}

#regions.contain.positions =
#	function(regions,positions)
#{
#	result=
#	sapply(1:nrow(regions),
#		function(r){
#			after.start = positions >= regions$start[r]
#			before.end = positions <= regions$end[r]
#			TRUE %in% (after.start & before.end)
#		}
#	)
#	return(result)
#}

regions.contain.positions =
	function(reg,positions)
{
	result=
	sapply(1:nrow(reg),
		function(r){
			after.start = positions >= reg$start[r]
			before.end = positions <= reg$end[r]
			as.character(reg$name[ after.start & before.end ])
		}
	)
	return(unlist(result))
}

both.positions.region =
	function(reg,p1p2)
{
	result=sapply(1:nrow(p1p2),
		function(i){
			regs.contain = reg$start <= p1p2[i,1] & reg$start <= p1p2[i,2] & reg$end >= p1p2[i,1] & reg$end >= p1p2[i,2]
			if(length(which(regs.contain))==0){regs.contain=NA}
			else{regs.contain=paste(reg$name[regs.contain],sep=',',collapse=',')}
			return(regs.contain)
		}
	)
	return(result)
}

empirical_distrib.regions =
	function(a,b,l,reg)
{
	if(a<1 | b<1 | l<1) return(0)
	# equal number of random positions
	rand.a = runif(a,1,l)
	rand.b = runif(b,1,l)
	# count number of regions that contain both a and b positions
	regions.a = regions.contain.positions(reg,rand.a)
	regions.b = regions.contain.positions(reg,rand.b)
	noverlap = length( intersect(regions.a,regions.b) )
	return(noverlap)
}

compute.MAST.MeDiChI.overlaps =
	function(masthits,chipfits,MeDiChI.pval.cutoff=0.05,overlap.distance=200,est.significance=FALSE,emp_iters=1000,comparison='distance',regions=NULL)
# navigate, translate, compare MAST hits table to MeDiChI data structure; count peak overlaps
{
	cat('comparing genome-wide MAST hits with p-value less than',MAST.pval.cutoff,'within',overlap.distance,'bp of MeDiChI peaks with p-value less than',MeDiChI.pval.cutoff,'\n')
	MAST.TFs=levels(factor(masthits$motif))
	overlaps=data.frame(seq=c(),mast.tf=c(),mast.pos=c(),chip.tf=c(),chip.exp=c(),chip.pos=c(),region=c())
	overlap.counts=data.frame(mast.tf=c(),chip.exp=c(),nmast=c(),nchip=c(),noverlap=c(),pval=c())

	for(mast.tf in MAST.TFs){
		cat('\n',mast.tf,'\n')
		for(chip.tf in names(chipfits)){
			for(exp in names(chipfits[[chip.tf]])){

				cat('\t',exp,'\n')
				seqs=names(chipfits[[chip.tf]][[exp]]$fits.fin)
				#'wg' stands for 'whole genome'
				nmast_wg=0
				nchip_wg=0
				noverlap_wg=0
				noverlap_random=list()

				# loop over each genome sequence contig that was fit separately
				for(chipseq in seqs){
					mastseq=NULL
					for(key in names(seqsyn)){
		#				cat('key',key,'\n')
						if(chipseq %in% seqsyn[[key]]){
							mastseq=key
							break
						}
					}
					if(is.null(mastseq)){
						cat('couldn\'t lookup key chipseq',chipseq,'\n')
						next
					} else {
						cat('\t\t',chipseq,'\n')
			#			cat('[MeDiChI]',chipseq,'[MAST]',mastseq,'\n')
					}

					noverlap=0
					seqmasthits=masthits[masthits$sequence==mastseq & grepl(mast.tf,masthits$motif),]
					hitpos=round((seqmasthits$hit_start+seqmasthits$hit_end)/2)
					nmast=length(hitpos)
		#			print(hitpos)

					sig.peaks=chipfits[[chip.tf]][[exp]]$fits.fin[[chipseq]][[1]]$coeffs.w.p.values
					# apply p-value filter to ChIP fit peaks
					sig.peaks=sig.peaks[sig.peaks[,3]<MeDiChI.pval.cutoff,]
					peak.positions=numeric()
					# special handling for 1-(or 0?) row results of filtering
					if(class(sig.peaks)=='numeric'){
						peak.positions=sig.peaks[1]
					}else{
						peak.positions=sig.peaks[,1]
					}
#					print(sig.peaks)
					npeaks=length(peak.positions)
					cat(nmast,'MAST hits vs.',npeaks,'MeDiChI peaks\n')

					if(nmast<1 | npeaks<1){
						noverlap=0
					}else if(comparison=='distance'){
						#noverlap=count_overlaps.distance(hitpos,peak.positions,dis=overlap.distance)
						distances=pairwise.site.distance(hitpos,peak.positions)
						overlapping=distances[distances$dis<overlap.distance,]
						if(nrow(overlapping)>0){
							# annotate overlaps contained in [promoter] regions
							ann=both.positions.region(regions[regions$seq==mastseq,],overlapping[,c(1,2)])
							overlaps=rbind(overlaps,data.frame(seq=mastseq,mast.tf=mast.tf,mast.pos=overlapping[,1],chip.tf=chip.tf,chip.exp=exp,chip.pos=round(overlapping[,2]),region=ann))
						}
						noverlap = nrow(overlapping)
					}else if(comparison=='regions'){
						# count number of regions that contain both MAST hits and ChIP peaks
						regions.hits = regions.contain.positions(regions[regions$seq==mastseq,],hitpos)
						regions.peaks = regions.contain.positions(regions[regions$seq==mastseq,],peak.positions)
						noverlap = length( intersect(regions.hits,regions.peaks) )
					}else{
						cat('invalid comparison type',comparison,'\n')
					}
		#			cat(' overlaps at dis',overlap.distance,'bp:',noverlap,'\n')
					nmast_wg = nmast_wg + nmast
					nchip_wg = nchip_wg + npeaks
					noverlap_wg = noverlap_wg + noverlap

					# estimation of statistical significance
					if(est.significance){
#						cat('estimating random empirical distribution of hit/peak overlaps...\n')
						if(comparison=='distance'){
							# if comparing actual hit/peak distances, then simulate random distribution (empirical)
							noverlap_random[[mastseq]]=
								replicate(emp_iters,empirical_distrib.distance(nmast,npeaks,seqlens[[mastseq]],dis=overlap.distance))
						}else if(comparison=='regions'){
							# count overlaps in regions for simulated random distributions of hits and peaks (emperical)
							noverlap_random[[mastseq]]=
								replicate(emp_iters,empirical_distrib.regions(nmast,npeaks,seqlens[[mastseq]],regions[regions$seq==mastseq,]))
						}else{
							cat('invalid comparison type',comparison,'\n')
						}
					}
				}

				# empirical p-value(?) fraction of time equal or greater overlaps in similar random distribution
				pval=NA
				if(est.significance){
					emp_overlaps=rowSums(do.call(cbind,noverlap_random))
					pval=length(which(emp_overlaps>=noverlap_wg))/length(emp_overlaps)
				}
				overlap.counts=rbind(overlap.counts,data.frame(mast.tf=mast.tf,chip.tf=chip.tf,chip.exp=exp,nmast=nmast_wg,nchip=nchip_wg,noverlap=noverlap_wg,pval=pval))
			}
		}
	}
	return(list(overlaps=overlaps,counts=overlap.counts))
}

#est.significance=TRUE
est.significance=FALSE
comparison='distance'
#comparison='regions'
emp_iters=250
#emp_iters=10
overlap.distance=100

overlaps=compute.MAST.MeDiChI.overlaps(masthits,MeDiChI.fits,MeDiChI.pval.cutoff=MeDiChI.pval.cutoff,overlap.distance=overlap.distance,est.significance=est.significance,emp_iters=emp_iters,comparison=comparison,regions=regions)

# add gene names
overlaps$overlaps = merge(overlaps$overlaps, coords[,c(1,2)],by.x='region',by.y='canonical_Name',sort=F,all.x=T)
# convert base chip TF names
overlaps$overlaps$chip.tf = sapply(overlaps$overlaps$chip.tf, conv.tf.name)

print(overlaps)


fname='MAST.MeDiChI.overlap.'
if(comparison=='distance'){
	fname=paste(fname,overlap.distance,'bp.',sep='')
}else if(comparison=='regions'){
	fname=paste(fname,'regions.',sep='')
}

write.table(overlaps$overlaps,paste(fname,'overlaps.tsv',sep=''),quote=F,row.names=F,sep='\t')
write.table(overlaps$counts,paste(fname,'overlap.counts.tsv',sep=''),quote=F,row.names=F,sep='\t')
