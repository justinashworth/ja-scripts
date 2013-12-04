#!/usr/bin/env python

# example usage:
# for f in [mast.txt's]; do [thisscript] $f > [somename]; done

import os,re,string,sys

class MASTGene:
	na_val = 1e10
	def  __init__(self,name='',pval=None,Eval=None):
		self.name = name
		if not pval: pval = self.na_val
		if not Eval: Eval = self.na_val
		self.pval = float(pval)
		self.Eval = float(Eval)
	def row(self,sep='\t'):
		return string.join([self.name, '%g'%self.pval, '%g'%self.Eval], sep)
	def __str__(self):
		return self.row()

''' Example result format in the mast.txt file:
VNG1097G

  LENGTH = 151  COMBINED P-VALUE = 2.36e-03  E-VALUE =      6.2
'''

re_score = re.compile('\n([^\n]+)\n\s+\n\s+LENGTH =\s+([0-9]+)\s+COMBINED P-VALUE\s+=\s+([0-9-+\.e]+)\s+E-VALUE =\s+([0-9-+\.e]+).+?DIAGRAM',re.DOTALL)

mastgenes = []
path = sys.argv[1]

for match in re_score.finditer(open(path).read()):
	name,length,pval,Eval = match.groups()
	mastgenes.append( MASTGene(name,pval,Eval) )
	#print str(mastgenes[len(mastgenes)-1])

print 'name\tpval\teval'
for mg in mastgenes:
	print mg.row()
