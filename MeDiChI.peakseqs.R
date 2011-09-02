# example for how to extract sequence regions around MeDiChI peaks
source('misc.R')

load('data/trh3.fits.RData')
# the MeDiChI Rdata file should result in a 'fits' object
peaks=MeDiChI.halo.peaks(fits,pval=0.05)

halo.genome=load.sequence('halo.genome.sequence.fa')

# get peak sequences
# upstream/downstream window for sequence around peak positions
window=c(-100,100)
seqs=peakseqs(peaks,halo.genome.seqs,window=window)
seqs=unlist(seqs)
# write to file
write.XStringSet(DNAStringSet(seqs),file=paste('Trh3.',window[1],'.',window[2],'.fa',sep=''))
