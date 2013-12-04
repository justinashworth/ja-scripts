#!/usr/bin/env python

import re,sys,string,os
from pssm import *
import cPickle

pssms = {}

for f in open(sys.argv[1]):
	mat = []
	for l in open(f.strip()):
		vals = l.strip().split()
		mat.append(vals)
#	print(mat)

	name = re.sub('.tsv','',f.strip())
	print name
	pssms[name] = pssm(biclusterName=name, pssm=mat, eValue='1e-100', nsites='100')

outFile = open('PredPSSMs.pkl','wb')
cPickle.dump(pssms,outFile)
outFile.close()

test = 'wt.pred.pal.ppm'
print(pssms[test].getConsensusMotif())
print(pssms[test].getMatrix())
print(pssms[test].getMemeFormatted())

test = '100.ALA.pal.ppm'
print(pssms[test].getConsensusMotif())
print(pssms[test].getMatrix())
print(pssms[test].getMemeFormatted())
