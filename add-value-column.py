#!/usr/bin/env python
import re,string,sys

from optparse import OptionParser
p=OptionParser()
p.add_option('-v','--values',help='values file ([id] [decimal#])')
p.add_option('-s','--sep',help='separator (default: \\t)',default='\t')
opt,args=p.parse_args()

#args specifies: destination file: [id] ..... 

values = {}
reval=re.compile('(.+)%s([0-9.-]+)'%opt.sep)
for match in re.finditer(reval,open(opt.values).read()):
	(id,val)=match.groups()
#	print id,val
	values[id]=val

dest=args[0]
outfile=open('%s.values'%dest,'w')
for line in open(dest):
	words=line.strip().split(opt.sep)
	id=words[0]
#	print id
	if not values.has_key(id):
		print 'warning, undefined id'
		continue
	val=values[id]
	lout = '%s%s%s\n' %(string.join([i for i in words[1:]],'\t'),opt.sep,val)
	outfile.write(lout)
outfile.close()
