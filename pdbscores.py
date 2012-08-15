#!/usr/bin/env python
from PDBFile import RosettaPDBFile
import os,sys,re,string

print 'file\tbound\tbinding'
pdblist = sys.argv[1]
for f in open(pdblist):
#	pdb = RosettaPDBFile( f, match_list = ['energy','pro_close','binding'] )
#	bound = pdb.E['energy'] - pdb.E['pro_close']
#	binding = pdb.E['binding']
#	print '%s\t%g\t%g' %(f, bound, binding)
	f = f.strip()

	pdb = RosettaPDBFile( f, match_list = ['energy','pro_close'] )
	bound = pdb.E['energy'] - pdb.E['pro_close']
	print '%s\t%g' %(f, bound)
