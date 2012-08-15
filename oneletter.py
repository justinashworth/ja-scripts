#!/usr/bin/env python

import os,re,sys,string
from AshworthUtil import threeletter
sep='\t'

fstr = open(sys.argv[1]).read()
for letter,three in threeletter.items():
	rethree = re.compile('%s%s' %(three,sep) )
	fstr = rethree.sub( '%s%s' %(letter,sep), fstr )
print fstr
