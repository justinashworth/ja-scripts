#!/usr/bin/env python

__author__ = "Justin Ashworth"
# jashworth@systemsbiology.org
# Institute for Systems Biology

import os,re,string,sys
from optparse import OptionParser
from Fasta import *

'''
DESCRIPTION
	given nucleotide sequences in fasta format, output all possible translated protein sequences in three frames for both directions
'''

class SixFrameTranslation(dict):
	def __init__(self,ntseq='',minlen=0):
		dict.__init__(self)
		self['fwd']=[None,None,None]
		self['rvs']=[None,None,None]
		self.translate(ntseq,minlen)

	def translate(self,ntseq,minlen=0):
		l = len(ntseq)
		for i in range(3):
			if i+3 >= l: continue
			self['fwd'][i] = translate_simple(ntseq[i:])
		rvs=rvs_comp_str(ntseq)
		for i in range(3):
			if i+3 >= l: continue
			self['rvs'][i] = translate_simple(rvs[i:])

	def __str__(self):
		out=[]
		for i in range(3):
			out.append( self['fwd'][i] )
		out.append('')
		for i in range(3):
			out.append( self['rvs'][i] )
		return string.join(out,'\n')

	def split(self,prefix,minlen):
		out=[]
		for dir in ['fwd','rvs']:
			for i in range(3):
				out.extend( protsplit(self[dir][i],'%s_%s'%(prefix,dir),minlen,i) )
		return out

def protsplit(protstr,prefix='',minlen=0,offset=0):
	re_prot = Fasta().re_protein
	out=[]
	for match in re_prot.finditer(protstr):
#		print prefix,offset,match.start(),match.end()
		if match.end() - match.start() < minlen: continue
		out.append( Fasta(protstr[match.start():match.end()],'%s_%i'%(prefix,match.start()*3+offset)) )
	return out

class App:
	def run(self,opt,args):
		seqs = FastaSeqs()
		seqs.go(args)
		for name,seq in sorted(seqs.seqs.items()):
			trans=SixFrameTranslation(seq)
#			print name
#			print str(trans)
#			print
			print string.join( [ s.fasta() for s in trans.split(name,opt.minlen) ], '\n' )

if __name__ == "__main__":
	p=OptionParser()
	p.add_option('--minlen',type='int',default=30)
	opt,args=p.parse_args()
	app = App()
	app.run(opt,args)
