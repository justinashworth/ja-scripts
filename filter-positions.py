#!/usr/bin/env python

'''
filter out position data that occurs within user-defined regions
'''

import os,re,string,sys
from optparse import OptionParser

class FilterPositions:
	
	def __init__(self):
		self.regions = {}

	def run(self,opt,args):
			self.read_regions(opt.regions,opt.debugout)
			self.filter_positions(args)

	def read_regions(self,file,debugout):
	
		# regions file format:
		#id	seq	strand	start	end
		#867	chr_1	-	300	10356
		#869	chr_1	-	17783	28449
		
		# stores regions by sequence key: script could be modified to avoid this layer
		sys.stderr.write('reading regions file %s\n' %file)
		for line in open(file):
			if line.startswith('id'): continue
			id,seq,strand,start,end = line.strip().split()
			if not self.regions.has_key(seq): self.regions[seq] = []
			self.regions[seq].append( (int(start),int(end)) )
		
		if debugout:
			for seq in self.regions:
				for start,end in self.regions[seq]:
					print '%s %i %i' %(seq,start,end)

	def filter_positions(self,files):

		# positions file format:
		#16150 12
		#16151 13
		#16152 15
		
		# looks for sequence key in filename
		re_chr = re.compile('(chr_[a-zA-Z0-9_]+)')
		for file in files:
			seq = re_chr.search(file).groups()[0]
			outfilename = '%s.filtered'%file
			sys.stderr.write('filtering positions file: %s (seq %s) outfile: %s\n' %(file,seq,outfilename))
			outfile = open(outfilename,'w')
			for line in open(file):
				pos,counts = line.strip().split()
				pos = int(pos)
				overlap = False
				for start,end in self.regions[seq]:
					if pos >= start and pos <= end:
						overlap = True
						break
				if not overlap: outfile.write(line)

if __name__ == "__main__":
	p=OptionParser()
	p.add_option('-r','--regions')
	p.add_option('--debugout',default=False,action='store_true')
	opt,args=p.parse_args()
	app = FilterPositions()
	app.run(opt,args)
