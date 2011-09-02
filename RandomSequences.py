#!/usr/bin/env python

__author__ = "Justin Ashworth"
# jashworth@systemsbiology.org
# Institute for Systems Biology

import os,re,string,sys,math,random
from optparse import OptionParser
from Fasta import Fasta

'''
DESCRIPTION
	Random sequence generator.
	Example:
		RandomSeqFactory -l 60 -n 1000 A 0.2 T 0.2 G 0.3 C 0.3
'''

class RandomSeqFactory:
	re_choice = re.compile('(\w+) +([0-9.]+)')

	def __init__(self,choicesstr=''):
		self.choices = {}
		self.read_choices_str(choicesstr)

	def read_choices_str(self,str):
		if str == '': return
		for m in self.re_choice.finditer(str):
			choice,prob = m.groups()
			prob = float(prob)
			if prob < 0 or prob > 1:
				print 'illegal prob %f' %prob; sys.exit()
			self.choices[choice] = prob
		sum = 0.0
		for p in self.choices.values():
			sum += p
		if abs(sum-1.0) >= 0.01:
			print 'probs do not sum to 1.0'; sys.exit()

	def read_choices_file(self,file):
		self.read_choices_str(open(file).read())

	def run(self,seqlen=10,nseqs=100,prefix=''):
		if self.choices == {}: return
		seqs = []
		for i in range(nseqs):
			seq = []
			for j in range(seqlen):
				cumprob = 0.0
				rnd = random.random()
				for choice,prob in self.choices.items():
					cumprob += prob
					if cumprob >= rnd:
						seq.append(choice)
						break
			seqs.append( Fasta( string.join(seq,''), prefix+str(i) ) )
		return seqs

if __name__ == "__main__":
	p=OptionParser()
	p.add_option('-f','--file')
	p.add_option('-l','--seqlen',type='int',default='10')
	p.add_option('-n','--nseqs',type='int',default='100')
	p.add_option('--prefix',default='random')
	opt,args=p.parse_args()
	app = RandomSeqFactory(string.join(args))
	if opt.file:
		app.read_choices_file(opt.file)
	seqs = []
	if opt.prefix:
		seqs = app.run(opt.seqlen,opt.nseqs,opt.prefix)
	else:
		seqs = app.run(opt.seqlen,opt.nseqs)
	for seq in seqs:
		print seq.table()
