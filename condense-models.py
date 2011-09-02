#!/usr/bin/env python
import os,re,string,sys
from optparse import OptionParser

'''
DESCRIPTION
	combines multiple exonic regions specified in a gff file for each given transcript into single continuous contigs
'''

#chr_1	JGI	exon	300	1153	.	-	.	name "fgenesh1_pg.C_chr_1000001"; transcriptId 867
re_exon=re.compile('(\w+)	JGI	exon	([0-9]+)	([0-9]+)	.	([+-])	.	name "[^"]+"; transcriptId ([0-9]+)')

class Model:
	sep='\t'
	def __init__(self,seq,start,end,strand,id):
		self.seq=seq
		self.start=int(start)
		self.end=int(end)
		self.strand=strand
		self.id=id
	def __str__(self):
		s='%s%s%i%s%i%s%s%s%s' %(self.seq,self.sep,self.start,self.sep,self.end,self.sep,self.strand,self.sep,self.id)
		return s
	def __cmp__(self,other):
		if self.seq < other.seq: return -1
		if self.seq > other.seq: return 1
		if self.start < other.start: return -1
		if self.start > other.start: return 1
		return -1

class App:
	def __init__(self):
		self.models={}
	def run(self,opt,args):
		for file in args:
			for m in re_exon.finditer(open(file).read()):
				seq,start,end,strand,id=m.groups()
		#		print seq,start,end,strand,id
				m=Model(seq,start,end,strand,id)
				if not models.has_key(id): models[id]=[]
				models[id].append(m)

		for id,ms in models.items():
			if len(ms) < 1: continue
			ms.sort()
		#	print id, ms[0].start, ms[0].end, ms[-1].start, ms[-1].end
			cmodel=ms[0]
			cmodel.end=ms[-1].end
			print cmodel

if __name__ == "__main__":
	p=OptionParser()
	opt,args=p.parse_args()
	app = App()
	app.run(opt,args)
