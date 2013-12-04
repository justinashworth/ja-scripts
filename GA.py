#!/usr/bin/env python

# TO DOs:
# support mixing of motifs with (slightly) different widths
#		may be necessary when using MEME motifs (MEME violates -maxw in -pal mode, for example)
# 	requires sliding/scoring for best alignment
#	generate center-gapped seed motifs
#		requires support for motifs of varying widths (above)
# cross-validation:
#		during training: at each step, leave 20% of targets out of training and base evolution on score for free set?
#		pre/post training: run whole training ten times with 10 random 10% leaveouts (merge ten final mots?)

# run with and without MEME seeds

import os, re, string, sys
from math import log
from DNAMotif import *
from Fasta import *
from optparse import OptionParser
from copy import deepcopy
import subprocess as sub
from scipy import stats
import random
import cPickle as pickle

FIMO_EXEC = 'fimo_4.9.0'
# total number of genes (for hypergeometric calcs)
N_ALL_GENES = 2400

################################################################################
# FUNCS

def write_motsfile(mots,fname=None):
	# write all mots to combined MEME file
	if not fname: fname = 'mots.MEME'
	outf = open(fname,'w')
	outf.write( mots[0].output_MEME_header() )
	for i in range(len(mots)):
		outf.write(mots[i].output_MEME())
	outf.close()

def write_motfile(mot,fname=None):
	# write all mots to combined MEME file
	if not fname: fname = 'mot.MEME'
	outf = open(fname,'w')
	outf.write( mot.output_MEME_header() )
	outf.write(mot.output_MEME())
	outf.close()

def run_FIMO(motpath,seqpath,bg=None,motname=None,thresh=0.0003):
	fimo_cmd = '%s -text --no-qvalue --norc --verbosity 1 --thresh %g' %(FIMO_EXEC,thresh)
	if bg: fimo_cmd = '%s -bgfile %s' %(fimo_cmd, bg)
	if motname: fimo_cmd = '%s -motif %s' %(fimo_cmd, motname)
	fimo_cmd = '%s %s %s' %(fimo_cmd, motpath, seqpath)
#	print fimo_cmd
	# capture output to stdout
	p = sub.Popen(fimo_cmd, stdout=sub.PIPE, stderr=sub.PIPE, shell=True)
	return p.communicate()[0]

def get_FIMO_genes(fimo):
	genes = []
	for l in fimo.split('\n'):
		if l.startswith('#'): continue
		gene = l.strip().split()
		if len(gene)<2: continue
		genes.append(gene[1])
	return set(genes)

def select_fimo_seq(fimo,targetgenes):
	# the first line is the header, and the last line is blank
	lines = fimo.split('\n')[1:-1]
	# selecting matched sequences from targetgenes only
	lines = [l for l in lines if l.split()[1] in targetgenes]
	if len(lines)<1: return None
	line = lines[ int(random.uniform(0,len(lines))) ]
	return line.strip().split()[7]

def phyper(genes,refgenes):
	inter = [g for g in genes if g in refgenes]
	tp = len(inter)
	trues = len(refgenes)
	samples = len(genes)
	return stats.hypergeom.sf(tp,N_ALL_GENES,trues,samples)

def score_overlap(genes,refgenes):
	inter = [g for g in genes if g in refgenes]
	tp = len(inter)
	trues = len(refgenes)
	samples = len(genes)
	# -log10 of phyper
#	score =  -1.0*log(stats.hypergeom.sf(tp,N_ALL_GENES,trues,samples),10)
	# adjust score based on number of refgenes matched (i.e. penalize over-fitting on a small fraction of refgenes)
#	score =  -1.0*log(stats.hypergeom.sf(tp,N_ALL_GENES,trues,samples),10) * float(tp)/trues
	# simple: precision * recall (roughly correlates with hypergeometric, but less sensitive to statistical noise)
	score = float(tp)/samples * (float(tp)/trues)
	return (score, tp, trues, samples)

def select_mot(mots):
	nmots = len(mots)
	inds = range(1,nmots+1)
	indcdf = [sum(inds[:i+1]) for i in range(len(inds))]
	rnd = int(round(random.uniform(0,sum(inds)-1)))+1
	topind = 0
	for ind in range(nmots):
#		print rnd, ind, indcdf[ind]
		if rnd >= indcdf[ind]: topind = ind
#	print 'selecting rank %i' %(topind+1)
	return mots[topind]

def center_gapped_motifs(mots,gap):
	gmots = []
	for mot in mots:
		for i in range(1,gap+1):
			gmots.append( center_gapped_motif(mot,i) )
	return gmots

def	pad_mots(mots,width):
	newmots = []
	for mot in mots:
		l = len(mot)
		if l >= width: continue
		pad = int(width)-l
		# even: add equal padding to both ends (in place)
		if pad%2==0:
			mot.pad_left(pad/2)
			mot.pad_right(pad/2)
			print 'padded mot %s: %s (%i)' %(mot.name,mot.consensus(),len(mot))
		# odd: pad half+1 right, and also return a new left-padded motif
		else:
			# (need to write a deep copy constructor for DNAMotif...)
			newmot = DNAMotif()
			newmot.set_width(width)
			newmot.matrix = mot.matrix
			newmot.data_type='probs'
			mot.pad_left((pad+1)/2)
			print 'padded mot %s: %s (%i)' %(mot.name,mot.consensus(),len(mot))
			newmot.pad_right((pad+1)/2)
			newmot.name = '%s_pad' %mot.name
			print 'padded mot %s: %s (%i)' %(newmot.name,newmot.consensus(),len(newmot))
			newmots.append(newmot)
		# (not implemented: other combinations of gaps and pads)
	for mot in newmots: mots.append(mot)

def find_MEME_motifs(seqfile,targetgenes,minw,maxw,nmotifs=2,bg=None,pal=False,maxseqs=100,prefix=''):
	# for MEME, limit to the 100 ChIP peaks with low p-values and highest intensities
	memegenes = []
	keep = min(len(targetgenes),maxseqs)
	counter = 0
	while len(memegenes)<keep:
		if counter >= len(targetgenes): break
		gene = targetgenes[counter][0]
		if not gene in memegenes: memegenes.append(gene)
		counter += 1
	print '%i genes for MEME' %len(memegenes)

	# write a fasta file for MEME targets
	seqs = FastaSeqs()
	seqs.loadseqs([seqfile])
	meme_fasta_path = '%smeme.fa' %prefix
	meme_fasta = open(meme_fasta_path,'w')
	for name,seq in seqs.seqs.items():
		if name in memegenes:
			meme_fasta.write(seq.fasta()+'\n')
	meme_fasta.close()

	# note that MEME will violate the -w and -maxw settings when -pal is used, potentially adding extra positions
	meme_cmd = 'meme_4.9.0 %s -bfile %s -minw %i -maxw %i -text -dna -nmotifs %i' %(meme_fasta_path,bg,minw,maxw,nmotifs)
	if pal: meme_cmd += ' -pal'
	print meme_cmd
	p = sub.Popen(meme_cmd, stdout=sub.PIPE, stderr=sub.PIPE, shell=True)
	meme_result = p.communicate()[0].split('\n')
	meme_mots = readMEME_from_lines(meme_result)
	for mot in meme_mots: mot.name = 'MEME_%s' %mot.name
	# force truncation of MEME motifs to maxw (-pal mode of MEME often violates its own -maxw parameter)
	for mot in meme_mots:
		l = len(mot)
		print mot.consensus(), l
		if l > maxw:
			mot.truncate(maxw)
			print 'MEME violated maxw: truncating to %s' %mot.consensus()
	return meme_mots

# END FUNCS
################################################################################
# SETUP

p=OptionParser()
p.add_option('-b','--bg')
p.add_option('-t','--targets', default='')
p.add_option('-c','--chip', default='')
p.add_option('-s','--seqs')
p.add_option('-g','--gap',type='int',default=0)
p.add_option('-m','--motifs_filelist')
p.add_option('-M','--MEME',action='store_true')
p.add_option('-x','--maxw',type='int',default=16)
p.add_option('-P','--pal',action='store_true')
p.add_option('-p','--pseudocounts',type='int',default='1')
p.add_option('--prefix',default='REFN')
p.add_option('-i','--iters',type='int',default=10)
opt,args=p.parse_args()

prefix = '%s/' %opt.prefix
if not opt.targets=='': prefix = '%s/%s/' %(prefix,opt.targets)
if not opt.chip=='': prefix = '%s/%s/' %(prefix,opt.chip)
os.system('mkdir -p %s' %prefix)

# load starting motifs
seed_mots = []
for f in open(opt.motifs_filelist):
#	print 'opening %s' %f.strip()
	mots = readMEME(f.strip())
	for mot in mots:
		print 'adding %s' %mot.name
		seed_mots.append(mot)
print 'loaded %i starting motifs' %len(seed_mots)

# weaken PWMs (sometimes appropriate for in silico or in vitro measurements with exaggerated IC?)
#for mot in seed_mots: mot.weaken(0.25)

if opt.gap > 0:
	print 'Adding center-gapped motifs'
	gapped_mots = center_gapped_motifs(seed_mots,opt.gap)
	print '%i center-gapped motifs' %len(gapped_mots)
	for mot in gapped_mots: seed_mots.append(mot)

targetgenes = []

if not opt.chip=='':
	chipgenes = []
	# load target genes (e.g. ChIP-bound genes)
	pval_cut = 0.1
	#position	intensity	pval	seq	expt	regionnames
	#932	1.4914345155365	0.791254125412541	Chr	trh3_2.2	NA
	for l in open(opt.chip):
		if l.startswith('position'): continue
		w = l.strip().split()
		gene = w[5]
		if gene=='NA': continue
		pval = float(w[2])
		if pval>pval_cut: continue
		intensity = float(w[1])
		chipgenes.append((gene,intensity))
	chipgenes.sort(key=lambda x: x[1])
	chipgenes.reverse()

	for gene,val in chipgenes:
		targetgenes.append(gene)
	targetgenes = set(targetgenes)

if not opt.targets=='':
	for t in open(opt.targets):
		targetgenes.append(t.strip())

print '%i target genes' %len(targetgenes)

# add MEME de novo motifs
if opt.MEME:
	minw=opt.maxw
	maxw=opt.maxw
	nmotifs=3
	pal=opt.pal
	maxseqs=100

	meme_mots = find_MEME_motifs(
		seqfile=opt.seqs,
		targetgenes=chipgenes,
		minw=minw,
		maxw=maxw,
		nmotifs=nmotifs,
		bg=opt.bg,
		pal=pal,
		maxseqs=maxseqs,
		prefix=prefix
	)
	print 'read %s MEME motifs' %len(meme_mots)
	for mot in meme_mots: seed_mots.append(mot)

# pad any short motifs
width = max(map(lambda x: len(x),seed_mots))
print 'padding all motifs to width %i' %width
pad_mots(seed_mots,width)

mots = deepcopy(seed_mots)
# TO DO diversify starting population of motifs

################################################################################
### begin refinement
MAX_MOTS = 100
CULL_FRAC = 0.4
TOP_CLONES = 1
NUM_CLONES = 1

# cache FIMO results (warning: potential memory issue!)
fimos = {}
# cache full mot scores and stats
scores = {}

topscore = None
nsame = 0

print 'Starting %s iters' %opt.iters
for iter in range(opt.iters):
	print '### iter %i ###' %iter

	if iter % 50 == 0:
		### save motifs
		pickle.dump(mots, open('%sPWMs.iter%i.p' %(prefix,iter),'wb'))

	any_new = False
	for mot in mots:
		if not mot.score: any_new = True

	if not any_new:
		print 'No new PWMs to score: stopping.'
		break

	nmots = len(mots)
	print 'Scoring %s motifs' %nmots

	# apparently FIMO wastes a bunch of time loading large MEME files, so this is not best
	#motpath = 'mots.MEME'
	#write_motsfile(mots,motpath)

	# 1. score all motifs for target enrichment
	# could be multiprocessed
	for mot in mots:
		# only score new/unscored motifs, since objective function is not changing
		if mot.score: continue
		motpath = '%smot.MEME' %prefix
		write_motfile(mot,motpath)
		# run FIMO for motif vs. all promoters
		fimo = run_FIMO(motpath,opt.seqs,opt.bg,mot.name)
		fimos[mot.name] = fimo
		genes = get_FIMO_genes(fimo)
		(score, tp, trues, samples) = score_overlap(genes, targetgenes)
		mot.score = score
		print '%s: %g' %(mot.name,mot.score)
		if not scores.has_key(mot.name): scores[mot.name] = (score, tp, trues, samples)

	mots.sort(key=lambda mot: mot.score)
	mots.reverse()

	# cache/write out full scores for whole population
	scorefile = open('%siter.%04i.scores' %(prefix,iter), 'w')
	for mot in mots:
		(score, tp, trues, samples) = scores[mot.name]
		scorefile.write('%s\t%g\t%i\t%i\t%i\n' %(mot.name,score,tp,trues,samples))
	scorefile.close()

	if iter==0:
		# save scores for seed_mots
		seed_mots = deepcopy(mots)

	report_leader = True
	if report_leader:
		leader = mots[0]
		# recompute score stats (not saved on a per-motif basis)
		(score, tp, trues, samples) = scores[leader.name]
		print 'Iter %i: Top motif %s %0.4g %i %i %i %s' %(iter, leader.name, score, tp, trues, samples, leader.consensus())
		if not topscore: topscore = score
		if score == topscore: nsame += 1
		if nsame >= 20:
			print 'Converged, stopping refinement'
			break
		if score > topscore:
			topscore = score
			nsame = 0

	# write out top motif at each iteration
	topmot = open('%srfn.top.iter.%04i.probs' %(prefix,iter), 'w')
	topmot.write( mots[0].output_default() )
	topmot.close()

	if iter+1 == opt.iters: break

	# 2. cull poorest-scoring motifs
	nkeep = int(MAX_MOTS*CULL_FRAC)
	nmots = len(mots)
	if nmots > nkeep:
		print 'culling %s low-scoring motifs' %(nmots-nkeep)
		del mots[nkeep:]

	nmots = len(mots)
	print '%s motifs survived' %nmots

	# 3. re-populate motif set

	# clone top motifs to enrich (increase likelihood of recombination with others)
	for i in range(TOP_CLONES):
		for j in range(NUM_CLONES):
			mots.append(mots[i])
	print 'cloned %ix%i top sequences' %(TOP_CLONES,NUM_CLONES)
	nmots = len(mots)
	# re-sort due to cloning via .append()
	mots.sort(key=lambda mot: mot.score)
	mots.reverse()

	# A.a every so often, spike with high-IC starting motifs in order to generate averaged/recombined motifs with current population
	max_mots = MAX_MOTS
	#if iter%5==0 and not iter==0 and not iter+1==opt.iters:
	if iter==10000:
		print 'iter %i: re-seeding with %i starting motifs' %(iter,len(seed_mots))
		for mot in seed_mots:
			mots.append(mot)
			# temporarily increase max_mots (since seeds are already scored)
		max_mots = max_mots + len(seed_mots)
		# re-sort due to seeding via .append()
		mots.sort(key=lambda mot: mot.score)
		mots.reverse()

	else:
		# A.b introduce matched sequences FROM TARGETS
		nemp = int( round( (max_mots-nmots) / 2 ))
		print 'adding %i empirically-updated motifs based on FIMO hits' %nemp
		for emp_i in range(nemp):
			# randomly select 5-20 FIMO match sequences (biased by motif rank)
			seqs = []
			nseqs = int(round(random.uniform(5,20)))
			while len(seqs) < nseqs:
				motname = select_mot(mots[:nmots]).name
				seq = select_fimo_seq(fimos[motname],targetgenes)
				if not seq: continue
				seqs.append(seq)
			empmot = DNAMotif(pseudocounts=opt.pseudocounts)
			empmot.probs_from_seqs(seqs)
			# mix with another existing motif
			empmot = mix_motifs(empmot, select_mot(mots[:nmots]))
			empmot.name = '%s_emp%i' %(iter,emp_i)
			mots.append(empmot)

	# B. add averaged child motifs
	navg = (max_mots-len(mots)) / 2
	print 'adding %i averaged motifs' %navg
	for child_i in range(navg):
		# pick parents, biased toward high rank
		mot1 = select_mot(mots[:nmots])
		mot2 = select_mot(mots[:nmots])
#		print str(mot1), str(mot2)
#		print 'mixing motifs %s and %s' %(mot1.name, mot2.name)
		# first half: average all parent positions
		child = mix_motifs(mot1,mot2,method='average')
		child.name = '%i_avg%i' %(iter,child_i)
#		# if a child is same as parent(s), set score as precomputed (it will not change)
		if child.same_probs(mot1):
			child.score = mot1.score
			child.name = mot1.name
		mots.append(child)

	# C. add recombined child motifs
	nrec = max_mots-len(mots)
	print 'adding %i recombined motifs' %nrec
	for child_i in range(nrec):
		mot1 = select_mot(mots[:nmots])
		mot2 = select_mot(mots[:nmots])
#		print 'mixing motifs %s and %s' %(mot1.name, mot2.name)
		# first half: average all parent positions
		child = mix_motifs(mot1,mot2,method='recombine')
		child.name = '%i_rec%i' %(iter,child_i)
		if child.same_probs(mot1):
			child.score = mot1.score
			child.name = mot1.name
		mots.append(child)

	# D. random mutations?

### output final motifs
pickle.dump(mots, open('%strained.PWMs.p' %prefix,'wb'))

# write top motif probfiles
details = open('%srfn.top.details' %prefix,'w')
N_PWMS = 10
top_mots = []
for mot in mots:
	if len(top_mots) == N_PWMS: break
	if mot.name in top_mots: continue
	top_mots.append(mot.name)
	fimo = fimos[mot.name]
	(score, tp, trues, samples) = scores[mot.name]
	p = phyper(get_FIMO_genes(fimo), targetgenes)
	details.write('%s %0.4g %i %i %i %g %s\n' %(mot.name, score, tp, trues, samples, p, mot.consensus()))
	pwmfile = open('%srfn.%s.probs' %(prefix,mot.name), 'w')
	pwmfile.write( mot.output_default() )
	pwmfile.close()
details.close()
#os.system('ls rfn*probs > %sprobfiles' %prefix)
#os.system('cd % && R --no-save < ../../motif.probfile.to.seqlogo.R' %prefix)
