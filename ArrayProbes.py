#!/usr/bin/env python

__author__ = "Justin Ashworth"
# jashworth@systemsbiology.org
# Institute for Systems Biology

import os,re,string,sys
from optparse import OptionParser

from Fasta import *
from AshworthUtil import rvs_comp_str

'''
This script creates probe sets for gene sequences
It can tile, or it can limit probes to specified regions (gene starts and stops, or exons with redundant ids) [--regions]
It can limit the number of probes per region supplied, starting from the 3' end [--probes-per-region] (for efficiently-distributed exon array)
It can also limit the number of probes per gene, starting from the 3' end [--probes-per-gene] (This is good for measuring either full-length transcripts, or poly-dT-amplified cDNA or cRNA)

*** TO DO ***:
add option to design anti-sense probes to expected transcript. This is critical if directly hybridizing total RNA to an array that isn't densely-tiled on both strands. Only methods involving opposite-strand synthesis, e.g. cRNA (Agilent Quick Amp) will hybridize correctly to a forward-strand array, (and vice versa: cRNA will not hybridize correctly to an antisense array)
'''

# custom comparison function for sorting probes
# could also have been a __cmp__ operator
re_probename=re.compile('([^_]+)_([^_]+)_([^_]+)_([^_]+)')
def probenamelt(p1,p2):
	m1=re_probename.match(p1)
	m2=re_probename.match(p2)
	if not m1 or not m2: return False
	seq1,start1,end1,str1=m1.groups()
	seq2,start2,end2,str2=m2.groups()
	if seq1 < seq2: return True
	if seq2 > seq1: return False
	if start1 < start2: return True
	if start1 > start2: return False
	if str1 < str2: return True
	return False

class Probe(Fasta):
	# a fasta sequence that is an array probe
	def __init__(self,start=0,strand='',parent='',**args):
		Fasta.__init__(self,**args)
		self.start=int(start)
		self.strand=strand
		self.parent=parent
	def end(self):
		return self.start+len(self)
	def __str__(self):
		if self.name != '':
			return '%s\t%s' %(self.name,self.seq.upper())
		else:
			return self.seq.upper()
	def __len__(self):
		return len(self.seq)
	# consider this 'static'
	def tableheader(self):
		header = ['name','parent','strand','start','end','sequence']
		return string.join(header,'\t')
	def table(self):
		row = [self.name,self.parent,self.strand,str(self.start),str(self.start+len(self)),self.seq.upper()]
		return string.join(row,'\t')
	def __lt__(self,other):
		# for sorting
		if probenamelt(self,other): return True
		if self.start < other.start: return True
		if self.start > other.start: return False
		if self.strand == '+' and other.strand == '-': return True
		return False

class Region:
	# sequence region (window)
	def __init__(self,seq,strand,start,end):
		self.seq = seq
		self.strand = strand
		self.start = int(start)
		self.end = int(end)
	def __cmp__(self,other):
		if self.seq < other.seq: return -1
		if self.seq > other.seq: return 1
		if self.seq != other.seq: return 0
		if self.strand == '+' and other.strand == '-': return -1
		if self.strand == '-' and other.strand == '+': return 1
		if self.strand != other.strand: return 0
		if self.end < other.end: return -1
		if self.end > other.end: return 1
		if self.start < other.start: return -1
		if self.start > other.start: return 1
		return 0
	def __str__(self):
		return '%s\t%s\t%i\t%i' %(self.seq,self.strand,self.start,self.end)

class RegionSet:
	# set of related Regions (e.g. exons for a single gene)
	def __init__(self,name=''):
		self.name = name
		self.regions = []
	def add(self,region):
		self.regions.append(region)
	def __str__(self):
		out = ['%s: RegionSet with %i Regions' %(self.name,len(self.regions))]
		for region in self.regions:
			out.append( str(region) )
		return string.join(out, '\n')

class ProbeSet:
	# collection of probes
	def __init__(self,probelength=60,gap=0,posoffset=0,negoffset=0,mask='',probespergene=0,probesperregion=0):
		self.probes=[]
		self.probelength=probelength
		self.gap=gap
		self.posoffset=posoffset
		self.negoffset=negoffset
		self.mask=mask
		self.probespergene=probespergene
		self.probesperregion=probesperregion
		self.regions={}

	def go(self,args):
		self.makeprobes(args)

	# limit probes to specific regions (e.g. genes, or exons)
	def read_regions(self,file):
		sys.stderr.write('reading regions from %s\n' %file)
		self.regions = {}
		# input format example
		#id seqname strand start end
		#867	chr_1	-	300	10356
		#867	chr_1	-	20000	21000
		#868	chr_2	-	23000	24000
		re_region = re.compile('([a-zA-Z0-9-_.]+)\s+([a-zA-Z0-9-_.]+)\s+([-+])\s+([0-9]+)\s+([0-9]+)')
		for line in open(file):
			m = re_region.match(line)
			if not m: continue
			id,seq,strand,start,end = m.groups()
			if strand.lower() in ['+','fwd']: strand = '+'
			elif strand.lower() in ['-','rvs']: strand = '-'
			else:
				print 'error unrecognized strand %s' %strand; sys.exit()
			start = int(start)
			end = int(end)
			# create a dict of RegionSets, keyed by id, for each (seq,strand)
			if not self.regions.has_key( (seq,strand) ): self.regions[(seq,strand)] = {}
			if not self.regions[(seq,strand)].has_key(id): self.regions[(seq,strand)][id] = RegionSet(id)
			self.regions[(seq,strand)][id].add( Region(seq,strand,start,end) )

		if self.regions == {}:
			# the regions might have lacked ids: try again with shorter regex
			# input format example
			#seqname strand start end
			#chr_1	-	300	10356
			#chr_1	-	20000	21000
			#chr_2	-	23000	24000
			re_region = re.compile('([a-zA-Z0-9-_.]+)\s+([-+])\s+([0-9]+)\s+([0-9]+)')
			for line in open(file):
				m = re_region.match(line)
				if not m: continue
				seq,strand,start,end = m.groups()
				#sys.stderr.write('%s %s %s %s' %(seq,strand,start,end))
				if strand.lower() in ['+','fwd']: strand = '+'
				elif strand.lower() in ['-','rvs']: strand = '-'
				else:
					print 'error unrecognized strand %s' %strand; sys.exit()
				id = '%s.%s.%s-%s' %(seq,strand,start,end)
				start = int(start)
				end = int(end)
				# create a dict of RegionSets, keyed by id, for each (seq,strand)
				if not self.regions.has_key( (seq,strand) ): self.regions[(seq,strand)] = {}
				if not self.regions[(seq,strand)].has_key(id): self.regions[(seq,strand)][id] = RegionSet(id)
				self.regions[(seq,strand)][id].add( Region(seq,strand,start,end) )

#		for seq,strand in self.regions:
#			print seq,strand
#			for regionset in self.regions[(seq,strand)].values():
#				print str(regionset)

	def makeprobes(self,fastafiles):
		f=FastaSeqs()
		f.loadseqs(fastafiles)
		print f.summarize()

		if self.regions == {}:
			for seq in f.seqs.values():
				# create region spanning whole seq
				l = len(seq)
				for d in ['+','-']:
					self.regions[(seq.name,d)] = {}
					ss = RegionSet('noid')
					ss.add( Region(seq.name,d,1,l) )
					self.regions[ (seq.name,d) ]['noid'] = ss
		else:
			# prevent lookup errors
			for seq in f.seqs.values():
				for d in ['+','-']:
					if not self.regions.has_key( (seq.name,d) ):
						self.regions[ (seq.name,d) ] = {'noid':RegionSet()}

		for seq in f.seqs.values():
			seqlen = len(seq)
			sys.stderr.write('%s\n' %seq.name)
			seqabbv=re.sub('_','',seq.name[ : min(len(seq.name),8) ])

			# the probes are created starting from the 3' end of the region. This makes it easy to create '3'-biased' probe sets of 'n' probes
			# forward strand
			for regionset in self.regions[(seq.name,'+')].values():
#				sys.stderr.write('+:\n%s\n' %str(regionset))
				# start with last region (exon), create probes from 3' end until nprobes reached
				nprobes_id = 0
				for region in reversed(sorted(regionset.regions)):
					nprobes_region = 0
					cursor = region.end - self.probelength + self.posoffset
					# note: strings 0-indexed (i), genome sequences 1-indexed (i+1)
					while cursor >= region.start:
						if cursor < 0 or cursor + self.probelength >= seqlen: break
						if self.probespergene > 0 and nprobes_id >= self.probespergene: break
						if self.probesperregion > 0 and nprobes_region >= self.probesperregion: break
						probeseq=seq.seq[cursor:cursor+self.probelength]
						if not re.search(opt.mask,probeseq): # skip probes containing masked sequence
							probename='%s_%i_+' %(seqabbv,cursor+1)
							p=Probe( start=cursor+1, strand='+', name=probename, seq=probeseq, parent=seq.name )
							self.probes.append(p)
							nprobes_id += 1
							nprobes_region += 1
						cursor -= self.probelength + self.gap

			# reverse strand
			for regionset in self.regions[(seq.name,'-')].values():
#				sys.stderr.write('-:\n%s\n' %str(regionset))
				nprobes_id = 0
				for region in sorted(regionset.regions):
					nprobes_region = 0
					cursor = region.start + self.negoffset
					while cursor < region.end - self.probelength:
						if cursor < 0 or cursor + self.probelength >= seqlen: break
						if self.probespergene > 0 and nprobes_id >= self.probespergene: break
						if self.probesperregion > 0 and nprobes_region >= self.probesperregion: break
						probeseq=seq.seq[cursor:cursor+self.probelength]
						if not re.search(opt.mask,probeseq):
							probename='%s_%i_-' %(seqabbv,cursor+1)
							p=Probe( start=cursor+1, strand='-', name=probename, seq=rvs_comp_str(probeseq), parent=seq.name )
							self.probes.append(p)
							nprobes_id += 1
							nprobes_region += 1
						cursor += self.probelength + self.gap

	def outputprobes(self,style='default'):
		out = ''
		if style == 'fasta':
			out = string.join([p.fasta() for p in self.probes],'\n')
		elif style == 'table':
			out = string.join([p.table() for p in self.probes],'\n')
		else: out = 'probe\tseq\n' + string.join([str(p) for p in self.probes],'\n')
		return out

class BlatMatch:
	# for parsing BLAT results
	#60	0	0	0	0	0	0	0	+	bd1x7_0_+	60	0	60	bd_1x7	26370	0	60	1	60,	0,	0,
	sep='\t'
	num='([0-9]+)' +sep
	fld='([^%s]+)' %sep +sep
	re_blat=re.compile( num*8 + '([+-])\t' + fld + num*3 + fld + num*3 )

	def __init__(self,match=None):
		self.match=0; self.mismatch=0; self.repmatch=0
		self.Ns=0
		self.Qgapcount=0; self.Qgapbases=0
		self.Tgapcount=0; self.Tgapbases=0
		self.strand=''
		self.Qname=''; self.Qsize=0; self.Qstart=0; self.Qend=0
		self.Tname=''; self.Tsize=0; self.Tstart=0; self.Tend=0
		# not used
#		self.blockcount=0, self.blockSizes=[], self.qStarts=[], self.tStarts=[]
		if match: self.parse(match)

	def parse(self,match):
		self.match,self.mismatch,self.repmatch,self.Ns,self.Qgapcount,self.Qgapbases,self.Tgapcount,self.Tgapbases,self.strand,self.Qname,self.Qsize,self.Qstart,self.Qend,self.Tname,self.Tsize,self.Tstart,self.Tend = match.groups()
		# int casting
		self.match=int(self.match); self.mismatch=int(self.mismatch); self.repmatch=int(self.repmatch)
		self.Ns=int(self.Ns)
		self.Qgapcount=int(self.Qgapcount); self.Qgapbases=int(self.Qgapbases)
		self.Tgapcount=int(self.Tgapcount); self.Tgapbases=int(self.Tgapbases)
		self.Qsize=int(self.Qsize); self.Qstart=int(self.Qstart); self.Qend=int(self.Qend)
		self.Tsize=int(self.Tsize); self.Tstart=int(self.Tstart); self.Tend=int(self.Tend)

	def perfect(self):
		# no partial
		if (self.match < self.Tsize) and (self.match < self.Qsize): return False
		# no mismatch
		if (self.mismatch != 0) or (self.repmatch !=0): return False
		# no gaps
		if (self.Qgapcount != 0) or (self.Tgapcount != 0): return False
		return True

	def good(self,mismatches=2,gaps=0):
		# no partial
		if (self.match < self.Tsize) and (self.match < self.Qsize): return False
		# mismatch limit
		if self.mismatch > mismatches: return False
		# gap limit
		if (self.Qgapbases > gaps) or (self.Tgapbases > gaps): return False
		return True

	def __str__(self):
		return string.join([str(x) for x in [self.match,self.mismatch,self.repmatch,self.Ns,self.Qgapcount,self.Qgapbases,self.Tgapcount,self.Tgapbases,self.strand,self.Qname,self.Qsize,self.Qstart,self.Qend,self.Tname,self.Tsize,self.Tstart,self.Tend] ], self.sep )

	def Rheader(self):
#		return string.join(['match','mismatch','repmatch','Ns','Qgapcount','Qgapbases','Tgapcount','Tgapbases','strand','Qname','Qsize','Qstart','Qend','Tname','Tsize','Tstart','Tend'],'\t')
		return string.join(['probe','sequence','strand','start','end'],'\t')
	def Routput(self):
		#return string.join([str(x) for x in [self.match,self.mismatch,self.repmatch,self.Ns,self.Qgapcount,self.Qgapbases,self.Tgapcount,self.Tgapbases,self.strand,self.Qname,self.Qsize,self.Qstart,self.Qend,self.Tname,self.Tsize,self.Tstart,self.Tend] ], '\t' )
		return string.join([str(x) for x in [self.Qname,self.Tname,self.strand,self.Tstart,self.Tend] ], '\t' )

if __name__ == "__main__":
	p=OptionParser()
	p.add_option('-l','--probelength',type='int',default=60)
	p.add_option('-g','--gap',type='int',default=0,help='can be negative for overlapping probes')
	p.add_option('--posoffset',type=int,default=0,help='offset for positive strands (from normal start position)')
	p.add_option('--negoffset',type=int,default=0,help='offset for negative strands')
	p.add_option('-o','--outfile')
	p.add_option('-f','--fastaout',default=False,action='store_true')
	p.add_option('-m','--mask',default='[acgtn]',help='regex pattern for masked sequence')
	p.add_option('-r','--regions',help='gene regions (or exons with redundant ids)')
	p.add_option('--probes-per-gene',type='int',default=0,help='maximum number of probes per gene id')
	p.add_option('--probes-per-region',type='int',default=0,help='maximum number of probes per region (e.g. exon)')
	opt,args=p.parse_args()

	app=ProbeSet(
		probelength=opt.probelength,
		gap=opt.gap,
		posoffset=opt.posoffset,
		negoffset=opt.negoffset,
		mask=opt.mask,
		probespergene=opt.probes_per_gene,
		probesperregion=opt.probes_per_region
	)

	if opt.regions: app.read_regions(opt.regions)
	app.go(args)

	style=''
	if opt.outfile:
		outf=open('%s'%opt.outfile,'w')
		outf.write( app.outputprobes() )
		outf.close()

		outf=open('%s.tsv'%opt.outfile,'w')
		outf.write( Probe().tableheader() + '\n' )
		outf.write( app.outputprobes('table') )
		outf.close()

		outf=open('%s.fa'%opt.outfile,'w')
		outf.write( app.outputprobes('fasta') )
		outf.close()
	else:
		if opt.fastaout: print app.outputprobes('fasta')
		else: print app.outputprobes()
