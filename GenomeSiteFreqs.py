#!/usr/bin/env python
import os,re,string,sys
from optparse import OptionParser
from Fasta import *
from AshworthUtil import all_combinations_gen,rvs_comp_str

'''
DESCRIPTION
	table of site frequencies from both strands in input sequences (fasta)
	such as for producing WEEDER background site Freq files
'''

class GenomeSiteFreqs:
	def __init__(self):
		self.sites = {}

	def initsites(self,length,pseudocount=1):
		sites = []
		choices = []
		for i in range(length):
			choices.append(['A','C','G','T'])
		all_combinations_gen(0,length,[],sites,choices)
		# pseudocounts (e.g. WEEDER)
		for site in sites: self.sites[site] = pseudocount

	def run_toplevel(self,opt,args):
		if not opt.sequences:
			print 'must supply sequence file (fasta)'; sys.exit()
		self.run(opt.sequences,opt.sitelength,opt.pseudocount)

	def run(self,fastafile,sitelength=6,pseudocount=1):
		sys.stderr.write('tabulating site frequencies for site length %i with pseudocount %i\n' %(sitelength,pseudocount))
		if not os.path.exists(fastafile):
			print 'sequence file %s not found' %fastafile; sys.exit()
		seqs = FastaSeqs()
		seqs.loadseqs([fastafile])
		self.initsites(sitelength,pseudocount)
		for seq in seqs.seqs.values():
			self.countsites(seq.seq,sitelength)

	def countsites(self,seq,length):
		n = range(len(seq)-length)
		for i in n:
			subseq = seq[i:i+length]
			self.sites[subseq] += 1
			self.sites[rvs_comp_str(subseq)] += 1

	def __str__(self):
		sites = self.sites.keys()
		sites.sort()
		out = []
		for site in sites:
			out.append( '%s %i' %(site,self.sites[site]) )
		return string.join(out, '\n')

if __name__ == "__main__":
	p = OptionParser()
	p.add_option('-s','--sequences')
	p.add_option('-l','--sitelength',type='int',default=6)
	p.add_option('-p','--pseudocount',type='int',default=1)
	opt,args = p.parse_args()
	app = GenomeSiteFreqs()
	app.run_toplevel(opt,args)
	print str(app)
