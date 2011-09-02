#!/usr/bin/env python
import os,re,string,sys
from optparse import OptionParser
p=OptionParser()
p.add_option('-d','--direction',default='')
p.add_option('-c','--lowcut',type='int',default=4)
opt,args=p.parse_args()

file = sys.argv[1]

#chr_1	1080	A	0	*	*	1	^F.	R	0	*	*	0	*	*	0	*	*	0	*	*	0	*	*
#chr_1	1081	T	0	*	*	1	.	[	0	*	*	0	*	*	0	*	*	0	*	*	0	*	*

for l in open(file):
	if not len(l) > 4: continue
	words = l.strip().split()
	seq = words[0]
	pos = words[1]
	l = len(words[3:])
	
	for i in range(l):
		
