#!/usr/bin/env python

'''
simple method to determine continuous regions in position data (such as RNA-seq mpileup reads-per-position)
assumes non-continuous position data, and that low signal positions have already been filtered out (faster to do this outside of python, e.g. grep)
combines across gaps (-g)
filters out small contigs (-m)
'''

import os,re,string,sys
from optparse import OptionParser

class Contigs:
	def run(self,opt,args):

		# input file format:
		#16150 12
		#16151 13
		#16172 15

		for file in args:
			positions = []
			for line in open(file):
				pos,counts = line.strip().split()
				positions.append(int(pos))

			positions = sorted(set(positions))

			contigs = []
			first = True
			start = 0
			lastp = 0
			for p in positions:
				if first:
					start = p
					first = False
				elif p > lastp+opt.gap+1:
					contigs.append( (start,lastp) )
					start = p
				lastp = p
			# final contig
			contigs.append( (start,lastp) )

			contigs = filter(lambda x: x[1]-x[0] >= opt.minlen, contigs)

			outfile = open('%s.contigs'%file,'w')
			for start,end in contigs:
				outfile.write('%s %s\n' %(start,end))

if __name__ == "__main__":
	p=OptionParser()
	p.add_option('-d','--debugout',default=False,action='store_true')
	p.add_option('-g','--gap',default='0',type='int')
	p.add_option('-m','--minlen',default='0',type='int')
	opt,args=p.parse_args()

	app = Contigs()
	app.run(opt,args)
