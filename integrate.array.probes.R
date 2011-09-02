# processed data from Bruz's pipeline
#ratios=read.delim('matrix_output.tsv')
ratios=read.delim('probe.ratios.tsv')

#ratios$probe=ratios$GENE
#ratios$GENE=NULL
#ratios$DESCRIPT=NULL

# strip out random probes
ratios=subset(ratios,!grepl('random',probe))
# strip out Agilent control probes(?)
ratios=subset(ratios,!grepl('^C',probe))

# add descriptive columns
ratios$strand='+'
ratios$strand[grep('-',ratios$probe)]='-'
# parse start positions from probe id
ratios$start=ratios$probe
ratios$start=as.numeric(gsub('.+_([0-9]+)_[+-]','\\1',ratios$start))
ratios$end[ratios$strand=='+'] = ratios$start[ratios$strand=='+']+60
# note: negative strand probes should also be annotated such that start<stop
#ratios$end[ratios$strand=='-'] = ratios$start[ratios$strand=='-']-60
ratios$end[ratios$strand=='-'] = ratios$start[ratios$strand=='-']+60
# parse sequence contig name from probe id
ratios$sequence=ratios$probe
ratios$sequence=gsub('(.+)_.+_[+-]','\\1',ratios$sequence)

# fix sequence names to correspond to reference sequence
ratios$sequence=gsub('chr','chr_',ratios$sequence)
ratios$sequence=gsub('chloropl','chloroplast',ratios$sequence)
ratios$sequence=gsub('mitochon','mitochondria',ratios$sequence)
ratios$sequence=gsub('chr_19a','chr_19a_19',ratios$sequence)
ratios$sequence=gsub('chr_19b','chr_19b_31',ratios$sequence)
ratios$sequence=gsub('chr_19c','chr_19c_29',ratios$sequence)

write_file=FALSE
if(write_file){
	write.table(ratios,'probe.values.tsv',quote=F,rownames=F,sep='\t')
}

#source('~/.Rgaggle')
#ds=getDatasetDescription()
#mat=subset(ratios,select=c('sequence','strand','start','end',"Day1.AM_v_Day3.PM.REF.ratio","Day1.PM_v_Day3.PM.REF.ratio","Day2.AM_v_Day3.PM.REF.ratio","Day2.PM_v_Day3.PM.REF.ratio","Day3.AM_v_Day3.PM.REF.ratio","Day3.PM_v_Day3.PM.REF.ratio","Day4.AM_v_Day3.PM.REF.ratio","Day4.PM_v_Day3.PM.REF.ratio","Day5.AM_v_Day3.PM.REF.ratio"))
#addTrack(ds,mat,type='quantitative.segment.matrix',name='Thaps.arrays')
