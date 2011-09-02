#!/usr/bin/env python
import os,re,string,sys
from optparse import OptionParser
from Fasta import *

'''
DESCRIPTION
	input: gff-formated gene model file
	output: fasta file containing reference sequences corresponding to models
'''

class GffToFasta:
	def run(self,opt,args):
		self.loadreference(opt.reference)
		self.findmodelseqs(opt.gff)
		self.writefile('%s.fa' %opt.gff)
	def loadreference(self,file):
		self.refseqs = FastaSeqs()
		self.refseqs.loadseqs([file])
	def findmodelseqs(self,file):
		self.modelseqs = []
		if not file or not os.path.exists(file):
			print 'invalid contigs file'; sys.exit()
		for l in open(file):
			fields = l.strip().split()
			if len(fields) < 7: continue
			seq = fields[0]
			start = int(fields[3])
			end = int(fields[4])
			direction = fields[6]
			if not self.refseqs.seqs.has_key(seq):
				sys.stderr.write('warning: seq %s is not in the reference sequence(s)\n' %seq)
				continue
			modelseq = self.refseqs.seqs[seq].subseq(start,end)
			if direction.lower() in ['rvs','-']: modelseq.rvscomp()
			modelseq.name = '%s %i-%i %s' %(seq,start,end,direction)
			self.modelseqs.append(modelseq)
	def writefile(self,filename):
		outfile = open(filename,'w')
		for modelseq in self.modelseqs:
			outfile.write(modelseq.fasta()+'\n')
		outfile.close()

if __name__ == "__main__":
	p=OptionParser()
	p.add_option('-r','--reference')
	p.add_option('-g','--gff')
	opt,args=p.parse_args()
	app = GffToFasta()
	app.run(opt,args)
