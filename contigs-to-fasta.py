#!/usr/bin/env python
import os,re,string,sys
from optparse import OptionParser
from Fasta import *

__author__ = 'Justin Ashworth'

'''
DESCRIPTION
	input: contigs (start, end) + fasta sequence(s)
	output: fasta file containing sequences corresponding to contigs
'''

class ContigsToFasta:

	def run(self,opt,args):
		self.loadSequences(opt.fasta)
		self.findContigSeqs(opt.contigs,opt.uppercase)

	def loadSequences(self,file):
		self.sequences = FastaSeqs()
		self.sequences.loadseqs([file])

	def findContigSeqs(self,file,uppercase):
		self.contigseqs = []
		if not file or not os.path.exists(file):
			print 'invalid contigs file'; sys.exit()
		#867	chr_1	-	300	10356
		#869	chr_1	-	17783	28449
		re_region = re.compile('([a-zA-Z0-9-_.]+)\s+([a-zA-Z0-9-_.]+)\s+([-+])\s+([0-9]+)\s+([0-9]+)')
		for m in re_region.finditer(open(file).read()):
			id,sequence,strand,start,end = m.groups()
			start = int(start)
			end = int(end)
			if not self.sequences.seqs.has_key(sequence):
				sys.stderr.write('warning: %s not found in the fasta sequence(s)\n' %sequence)
				continue
			contigseq = self.sequences.seqs[sequence].subseq(start,end)
			if uppercase: contigseq.seq = contigseq.seq.upper()
			if strand == '-': contigseq.rvscomp()
			#contigseq.name = '%s_%i-%i_%s' %(sequence,start,end,strand)
			contigseq.name = '%s:%s_%i-%i_%s' %(id,sequence,start,end,strand)
			self.contigseqs.append(contigseq)

	def __str__(self):
		out = []
		for contig in self.contigseqs:
			out.append(contig.fasta())
		return string.join(out,'\n')

if __name__ == "__main__":
	p=OptionParser()
	p.add_option('-f','--fasta')
	p.add_option('-c','--contigs')
	p.add_option('-U','--uppercase',default=False,action='store_true')
	opt,args=p.parse_args()
	app = ContigsToFasta()
	app.run(opt,args)
	print str(app)
