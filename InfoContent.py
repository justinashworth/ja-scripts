#!/usr/bin/env python

__author__ = "Justin Ashworth"
# jashworth@systemsbiology.org
# Institute for Systems Biology

import os,re,string,sys,math
from optparse import OptionParser

class InfoContent:
	def __init__(self,base=2,background_entropy=2):
		self.base=base
		self.bge=background_entropy
	def run(self,seqs):
		maxlen = max([len(seq) for seq in seqs])
		ics = []
		for i in range(maxlen):
			column = []
			for seq in seqs:
				if i >= len(seq):	continue
				column.append(seq[i])
			ics.append(self.ic(column))
		return ics
	def ic(self,symbols):
		counts = {}
		for s in symbols:
			if not counts.has_key(s): counts[s]=0
			counts[s] += 1
		ic = self.bge
		l = len(symbols)
		for s,c in counts.items():
			prob = float(c)/l
			ic += prob * math.log(prob,self.base)
		return ic

if __name__ == "__main__":
	p=OptionParser()
	opt,args=p.parse_args()
	app = InfoContent()
