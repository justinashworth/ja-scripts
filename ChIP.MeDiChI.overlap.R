# cross-compare peak overlaps for MeDiChI fits on FY's new Lrp ChIP data

# for pairwise distance funcs:
source('peak.overlap.R')

# bootstrap p-value for MeDiChI fits
MeDiChI.pval.cutoff=0.05
#MeDiChI.pval.cutoff=0.02

seqlens=list(
	'HALCHR'=2014239,
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

# regions in which to count overlaps
if(!exists('regions')){
	upstream=400
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

annotate.regions =
	# SLOOOOW (because R sucks at for loops)
	function(positions,regions)
{
	result=
	sapply(positions,
		function(x){
			for(i in 1:nrow(regions)){
				if(x>regions$start[i] & x<regions$end[i]){
					return(as.character(regions$name[i]))
				}
			}
			return(NA)
		}
	)
	return(unlist(result))
}

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

regions.contain.positions =
	function(regions,positions)
{
	result=
	sapply(1:nrow(regions),
		function(r){
			after.start = positions >= regions$start[r]
			before.end = positions <= regions$end[r]
			TRUE %in% (after.start & before.end)
		}
	)
	return(result)
}

empirical_distrib.regions =
	function(a,b,l,regions)
{
	if(a<1 | b<1 | l<1) return(0)
	# equal number of random positions
	rand.a = runif(a,1,l)
	rand.b = runif(b,1,l)
	# count number of regions that contain both a and b positions
	regions.a = which( regions.contain.positions(regions,rand.a) )
	regions.b = which( regions.contain.positions(regions,rand.b) )
	noverlap = length( intersect(regions.a,regions.b) )
	return(noverlap)
}

MeDiChI.peakpos =
	function(peaktab,pval=0.05)
{
	# apply p-value filter to ChIP fit peaks
	peaks=peaktab[peaktab[,3]<MeDiChI.pval.cutoff,]
	peak.positions=numeric()
	# special handling for 1-(or 0?) row results of filtering
	if(class(peaks)=='numeric'){
		peak.positions=peaks[1]
	}else{
		peak.positions=peaks[,1]
	}
	return(peak.positions)
}

# specific to latest data combined MeDiChI data structure ( {tf}{experiment}$fits.fin... )
# first list level is unnecessary (second list level names (experiment) are unique)
MeDiChI.exps = unlist(MeDiChI.fits,recursive=FALSE)

compute.MeDiChI.overlaps =
	function(chipfits,MeDiChI.pval.cutoff=0.05,overlap.distance=200,est.significance=FALSE,emp_iters=1000,comparison='distance',regions=NULL)
# navigate, translate, compare MAST hits table to MeDiChI data structure; count peak overlaps
{
	cat('comparing genome-wide overlaps within',overlap.distance,'bp for MeDiChI peaks with p-value less than',MeDiChI.pval.cutoff,'\n')
	overlaps=data.frame(chip1=c(),chip2=c(),n1=c(),n2=c(),noverlap=c(),pval=c())

	for(i in 1:length(chipfits)){
		exp1 = names(chipfits)[i]
		for(j in i:length(chipfits)){
			exp2 = names(chipfits)[j]
			cat('\n',exp1,exp2,'\n')

			#'wg' stands for 'whole genome'
			n1_wg=0
			n2_wg=0
			noverlap_wg=0
			noverlap_random=list()

			for(seq in names(seqlens)){
				cat(seq,'\n')
				noverlap=0

#				if( !seq %in% chipfits[[exp1]]$fits.fin | !seq %in% chipfits[[exp1]]$fits.fin
				peaks1=MeDiChI.peakpos( chipfits[[exp1]]$fits.fin[[seq]][[1]]$coeffs.w.p.values, MeDiChI.pval.cutoff )
				peaks2=MeDiChI.peakpos( chipfits[[exp2]]$fits.fin[[seq]][[1]]$coeffs.w.p.values, MeDiChI.pval.cutoff )

				n1=length(peaks1)
				n2=length(peaks2)
				cat(seq,n1,n2,'\n')

				if(n1<1 | n2<1){
					noverlap=0
				}else if(comparison=='distance'){
					noverlap=count_overlaps.distance(peaks1,peaks2,dis=overlap.distance)
				}else if(comparison=='regions'){
					# count number of regions that contain both MAST hits and ChIP peaks
					regions1 = which( regions.contain.positions(regions,peaks1) )
					regions2 = which( regions.contain.positions(regions,peaks2) )
					noverlap = length( intersect(regions1,regions2) )
				}else{
					cat('invalid comparison type',comparison,'\n')
				}
				#cat(' overlaps at dis',overlap.distance,'bp:',noverlap,'\n')
				n1_wg = n1_wg + n1
				n2_wg = n2_wg + n2
				noverlap_wg = noverlap_wg + noverlap

				# estimation of statistical significance
				if(est.significance){
					#cat('estimating random empirical distribution of hit/peak overlaps...\n')
					if(comparison=='distance'){
						# if comparing actual hit/peak distances, then simulate random distribution (empirical)
						noverlap_random[[seq]]=
							replicate(emp_iters,empirical_distrib.distance(n1,n2,seqlens[[seq]],dis=overlap.distance))
					}else if(comparison=='regions'){
						# count overlaps in regions for simulated random distributions of hits and peaks (emperical)
						noverlap_random[[seq]]=
							replicate(emp_iters,empirical_distrib.regions(n1,n2,seqlens[[seq]],regions))
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
			overlaps=rbind(overlaps,data.frame(chip1=exp1,chip2=exp2,n1=n1_wg,n2=n2_wg,noverlap=noverlap_wg,pval=pval))
			# flip chip ids to fill symmetric matrix (useful later for xtabs, heatmaps, etc)
			# (except for self-self)
			if(!chip1==chip2){
				overlaps=rbind(overlaps,data.frame(chip1=exp2,chip2=exp1,n1=n1_wg,n2=n2_wg,noverlap=noverlap_wg,pval=pval))
			}
		}
	}
	return(overlaps)
}

est.significance=TRUE
#est.significance=FALSE
comparison='distance'
#comparison='regions'
emp_iters=1000
#emp_iters=250
#emp_iters=10
overlap.distance=200

overlaps=compute.MeDiChI.overlaps(MeDiChI.exps,MeDiChI.pval.cutoff=MeDiChI.pval.cutoff,overlap.distance=overlap.distance,est.significance=est.significance,emp_iters=emp_iters,comparison=comparison,regions=regions)
print(overlaps)

fname='MeDiChI.overlap.'
if(comparison=='distance'){
	fname=paste(fname,overlap.distance,'bp.',sep='')
}else if(comparison=='regions'){
	fname=paste(fname,'regions.',sep='')
}
fname=paste(fname,'tsv',sep='')

write.table(overlaps,fname,quote=F,row.names=F,sep='\t')
