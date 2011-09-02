#!/usr/bin/env python
__author__ = "Justin Ashworth"
import os,re,string,sys
from optparse import OptionParser
from AshworthUtil import ChainPos

'''
DESCRIPTION
Work in progress. Basic alignment representation to support simple tasks.

Example of alignment file:

template P.OT3.FL11    A.22,APLREISKITGLAESTIHERIRKLRESG B.22,APLREISKITGLAESTIHERIRKLRESG
model    Halo.AsnC          KPFSEIAREIDMSSATVHDRVSRLEEAG      KPFSEIAREIDMSSATVHDRVSRLEEAG
model    Halo.Trh2          TSAPDIADEVNVSAGTIRNRITHLEAAG      TSAPDIADEVNVSAGTIRNRITHLEAAG
'''

class ChainSeq:
	def __init__(self,start,sequence):
		self.start = start
		self.sequence = sequence

	def __str__(self):
		return '%s,%s' %(self.start,self.sequence)

class AlignmentModel:
	def __init__(self):
		self.type = None
		self.name = None
		self.chainseqs = []

	def parseline(self,line):
		words=line.strip().split()
		if len(words) < 3:
			print 'error parsing AlignmentModel from line %s' %line; return
		self.type = words[0]
		if self.type != 'model':
			print 'error, not tagged as model: %s' %line; return
		self.name = words[1]
		for seq in words[2:]:
			self.chainseqs.append(seq)

	def __str__(self):
		out = '%10s %16s' %(self.type,self.name)
		for cs in self.chainseqs: out += '%6s' %'' + cs
		return out

class AlignmentTemplate:
	def __init__(self):
		self.type = None
		self.name = None
		self.chainseqs = []

	def parseline(self,line):
		words=line.strip().split()
		if len(words) < 3:
			print 'error parsing AlignmentTemplate from line %s' %line; return
		self.type = words[0]
		if self.type != 'template':
			print 'error, not tagged as template: %s' %line; return
		self.name = words[1]
		for chainseq in words[2:]:
			start,seq = chainseq.split(',')
			chain,pos = start.split('.')
			start = ChainPos(chain,pos)
			self.chainseqs.append( ChainSeq(start,seq) )

	def __str__(self):
		out = '%10s %16s' %(self.type,self.name)
		for cs in self.chainseqs: out += ' ' + str(cs)
		return out

	# temporary hack that might have better OO inheritance solution?
	def model(self):
		m = AlignmentModel()
		m.type = self.type
		m.name = self.name
		for cs in self.chainseqs:
			m.chainseqs.append( cs.sequence )
		return m

class Alignment:
	def __init__(self):
		self.templates = {}
		self.models = {}

	def run_toplevel(self,opt,args):
		if not opt.alignfile:
			print 'error, must specify alignfile (-f option)'; return
		self.run(opt.alignfile)

	def run(self,alignmentfile):
		self.loadmodels(alignmentfile)

	def loadmodels(self,alignfile):
		if not alignfile or not os.path.exists(alignfile):
			print 'error, %s doesn\'t exist' %alignfile; return
		for line in open(alignfile):
			type = line.strip().split()[0]
			if type == 'template':
				template = AlignmentTemplate()
				template.parseline(line)
				self.templates[template.name] = template
			elif type == 'model':
				model = AlignmentModel()
				model.parseline(line)
				self.models[model.name] = model

	def __str__(self):
		out = []
		for t in self.templates.values():
			out.append(str(t))
		for m in self.models.values():
			out.append(str(m))
		return string.join(out,'\n')

if __name__ == "__main__":
	p=OptionParser()
	p.add_option('-f','--alignfile')
	opt,args=p.parse_args()
	app = Alignment()
	app.run(opt,args)
	print str(app)
