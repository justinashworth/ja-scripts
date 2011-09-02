#!/usr/bin/env python
import os,re,string,sys
from optparse import OptionParser
p=OptionParser()
p.add_option('-R','--Routput',default=False,action='store_true')
opt,args=p.parse_args()

from ArrayProbes import *

#allmatches={}
if opt.Routput: print BlatMatch().Rheader()
for file in args:
	for line in open(file):
		match = BlatMatch().re_blat.match(line)
		# pass through non-matching (header) lines
		if not match:
			if not opt.Routput:
				print line,
			continue
		bm = BlatMatch(match)
#		if not bm.good(): continue
		if not bm.perfect(): continue
		# just print. otherwise too much memory for such a simple task
		if opt.Routput: print bm.Routput()
		else: print str(bm)
#		if not allmatches.has_key(bm.Qname): allmatches[bm.Qname] = []
#		allmatches[bm.Qname].append(bm)
	
#filteredmatches={}
#for probe,matches in allmatches.items():
#	print probe
#	perfect=[]
	# these T/F flags specific to distinguishing between assembled and "bottom drawer" contigs for Thaps
#	chr_match=False
#	bd_match=False
#	org_match=False
#	for match in matches:
#		print str(match)
#		if not match.perfect(): continue
#		if not match.good(): continue
#		if not filteredmatches.has_key(probe): filteredmatches[probe]=[]
#		filteredmatches[probe].append(match)
#		print probe,'perfect match'
		# Tname is the name of genome contig
#		if re.match('bd_',match.Tname): bd_match=True
#		if re.match('chr_',match.Tname): chr_match=True
#		if re.match('chloroplast',match.Tname): org_match=True
#		if re.match('mitochondria',match.Tname): org_match=True
#		perfect.append(match)
#	keepers=[]
#	if bd_match and (chr_match or org_match):
#		for match in perfect:
#			if re.match('bd_',match.Tname): continue
#			keepers.append(match)
#	else: keepers = perfect
#	filteredmatches[probe]=keepers

#for matches in filteredmatches.values():
#for matches in allmatches.values():
#	print string.join([str(m) for m in matches], '\n')
