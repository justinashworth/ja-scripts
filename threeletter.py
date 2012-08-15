#!/usr/bin/env python

import os,re,sys,string
from AshworthUtil import threeletter
sep='\t'

fstr = open(sys.argv[1]).read()
for letter,three in threeletter.items():
	relt = re.compile('%s%s%s' %(sep,letter,sep) )
	fstr = relt.sub( '%s%s%s' %(sep,three,sep), fstr )
print fstr
