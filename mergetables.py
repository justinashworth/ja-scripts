#!/usr/bin/env python
import os,re,sys,string

from optparse import OptionParser
p=OptionParser()
p.add_option('-s','--sep',help='separator (default: \\t)',default='\t')
p.add_option('-o','--out',help='outfile name',default='mergedfile')
p.add_option('-c','--chrseg',help='use match chromosome segment definition as key',action='store_true')
opt,args=p.parse_args()

class ChrSeg:
	rematch=re.compile('([Cc]hr[0-9]+[a-z]*)\s+([\+.-])\s+([0-9.-]+)\s+([0-9.-]+)')
	def __init__(self):
		self.chr = ''
		self.strand = '.'
		self.start = 0
		self.end = 0
	def parse(self,match):
		(chr,strand,start,end)=match.groups()
		self.chr=chr
		self.strand=strand
		self.start=int(start)
		self.end=int(end)
	def __str__(self):
		return '%s%s%s%s%i%s%i' %(self.chr,opt.sep,self.strand,opt.sep,self.start,opt.sep,self.end)
	def __eq__(self,other):
		return self.chr == other.chr and self.strand == other.strand and self.start == other.start and self.end == other.end
	def __hash__(self):
		return hash(str(self))

# combine data columns from mutiple files that share a common first column if ids, that may be out of order or non-overlapping

values={}
# lists to preserve key order for output from dicts
roworder=[]
columnkeys=[]

if len(args)==0: sys.exit()
for file in args:
	if not os.path.exists(file):
		print 'skipping %s' %file
		continue
	# first column is primary key
	for line in open(file):
		words=line.strip().split(opt.sep)
		id=words[0]
		# special case (multicolumn key, see class above)
		if opt.chrseg:
			m=ChrSeg.rematch.match(line)
			if m:
				id=ChrSeg()
				id.parse(m)
			else:
				print 'warning, failed ChrSeg match'
				continue

		if not values.has_key(id):
			values[id]={}
			roworder.append(id) # preserve original key order as encountered
		starti=1
		if opt.chrseg: starti=4 # (special case)
		for i in range(len(words[starti:])):
			key='%s%s'%(file,str(i)) # secondary key: file-specific column index
			if not key in columnkeys: columnkeys.append(key) # preserve key (column) order
			values[id][key]=words[i+starti]

outf=open(opt.out,'w')
for id in roworder:
	line=[str(id)]
	for key in columnkeys:
		if not values[id].has_key(key): line.append('') # prevent column collapse for missing values
		else: line.append( values[id][key] )
	outf.write('%s\n'%string.join(line,opt.sep))
outf.close()
