#!/usr/bin/env python

import os, re,string,sys

names = sys.argv[1:]

by = 1
splitchr = '.'
field = 0

ops = []

re_num = re.compile('([0-9]+)')

for n in names:
	match = re_num.search(n)
	new = n
	if match:
		ind = int(match.groups()[0])
		newind = ind+by
		new = re_num.sub('%04i'%newind,n)
	ops.append('mv %s renum.%s' %(n,new))

for op in ops:
	os.system(op)
