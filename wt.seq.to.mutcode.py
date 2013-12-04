#!/usr/bin/env python

import os,sys,re,string

seq = open(sys.argv[1]).read()

count = 1

print 'mutation	pos	wt	mut'
for letter in seq.strip():
	print '%s%i%s\t%i\t%s\t%s' %(letter,count,letter,count,letter,letter)
	print '%s%i%s\t%i\t%s\t%s' %(letter,count,'X',count,letter,'X')
	count += 1
