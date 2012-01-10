#!/usr/bin/env python
import os,re,string,sys
from optparse import OptionParser

'''
DESCRIPTION
	Doesn't actually map sequences--coordinates for probe->genome mapping should be known already (e.g. using BLAT).
	This just handles the numerical bookkeeping of probe(sequence,strand,start,end)->gene(sequence,strand,start,end)
	file format example for probes and genes:
	probe/gene	seq	strand	start	end
	chr11a_538659_+	chr_11a	+	538658	538718
	chr11a_538799_+	chr_11a	+	538798	538858
'''

class Probe:
	#chr11a_538659_+ chr_11a + 538658  538718
	re_probe = re.compile('([A-Za-z0-9_+-.]+)\t([^\t]+)\t([+-])\t([0-9]+)\t([0-9]+)')
	def __init__(self,name,seq,strand,start,end):
		self.name = name
		self.seq = seq
		self.strand = strand
		self.start = start
		self.end = end
	def __str__(self):
		out = [self.name,self.seq,self.strand,str(self.start),str(self.end)]
		return string.join(out,'\t')
	def __cmp__(self,other):
		# sorting operator
		if self.seq < other.seq: return -1
		if self.seq > other.seq: return 1
		if self.strand == '+' and other.seq =='-': return -1
		if self.strand == '-' and other.seq =='+': return 1
		if self.start < other.start: return -1
		if self.start > other.start: return 1
		if self.end < other.end: return -1
		if self.end > other.end: return 1
		return 0

class Gene:
	#20563 chr_1 - 33946 34494
	#31076 chr_1 + 34958 36758
	re_gene = re.compile('([A-Za-z0-9_+-.]+)\t([^\t]+)\t([+-])\t([0-9]+)\t([0-9]+)')
	def __init__(self,id,seq,strand,start,end):
		self.id = id
		self.seq = seq
		self.strand = strand
		self.start = start
		self.end = end
		self.probes = []
	def setProbes(self,probes=[]):
		self.probes=probes
	def __str__(self):
		out = [self.id,self.seq,self.strand,str(self.start),str(self.end)]
		out = string.join(out,'\t')
		out += ': %s probes' %len(self.probes)
		return out
	def probemap(self):
		for probe in self.probes:
			return '%s\t%s' %(self.probe,self.id)

class App:
	def __init__(self):
		self.genes = {}

	def run(self,opt,args):
		self.loadGenes(opt.genes)
		sys.stderr.write(self.geneOverlap()+'\n')
		self.loadProbes(opt.probes,opt.antisense)
		if opt.threeprime:
			self.threePrime(opt.threeprime)

	def loadGenes(self,genes):
		print 'probe\tsequence\tstrand\tstart\tend\tgeneid'
		self.genes = {}
		re_gene = Gene.re_gene
		gcount = 0
		for match in re_gene.finditer(open(genes).read()):
			id,seq,strand,start,end = match.groups()
			start = int(start)
			end = int(end)
			if not self.genes.has_key(seq): self.genes[seq] = {'+':[],'-':[]}
			self.genes[seq][strand].append( Gene(id,seq,strand,start,end) )
			gcount += 1
		sys.stderr.write('%s genes read\n' %gcount)

	def geneOverlap(self):
		# compare gene boundaries to identify overlaps
		# (some sets of annotations contain redundancies)
		out = []
		for seq in self.genes:
#			print seq
			for strand in self.genes[seq]:
#				print strand
				genes=self.genes[seq][strand]
				l=len(genes)
				for i in range(l):
					gene1=genes[i]
					for j in range(i+1,l):
						gene2=genes[j]
						if (gene2.start >= gene1.start and gene2.start <= gene1.end) \
						  or (gene2.end >= gene1.start and gene2.end <= gene1.end):
							if gene2.start >= gene1.start and gene2.end <= gene1.end:
								out.append('gene %s is contained within gene %s' %(gene2.id,gene1.id))
							elif gene1.start >= gene2.start and gene1.end <= gene2.end:
								out.append('gene %s is contained within gene %s' %(gene1.id,gene2.id))
							else:
								out.append('gene %s overlaps gene %s' %(gene1.id,gene2.id))
		return string.join(out,'\n')

	def loadProbes(self,probes,antisense=False):
		re_probe = Probe.re_probe
		pcount = 0
		pcount_assigned = 0
		pcount_genes = 0
		for match in re_probe.finditer(open(opt.probes).read()):
			name,seq,strand,start,end = match.groups()
			pcount += 1
			if antisense:
				if strand == '+': strand = '-'
				elif strand == '-': strand = '+'
				else: print 'unknown strand!'
			start=int(start)
			end=int(end)
			# correct for flipped start and end coords
			if start > end:
				sys.stderr.write('correcting flipped start/end\n')
				newstart = end
				end = start
				start = newstart
			# assign probe to overlapping gene region(s)
			halfprobe = abs(end-start)/2
			if not self.genes.has_key(seq):
				sys.stderr.write('sequence \'%s\' not in genes!\n' %seq)
				continue
			for gene in self.genes[seq][strand]:
#				print str(gene)
				# both probes and genes should be numbered with stop>start, regardless of strand
#				print start,end,gene.start,gene.end
				if start + halfprobe - 1 >= gene.start and end <= gene.end + halfprobe - 1:
					pcount_assigned += 1
					if len(gene.probes) < 1: pcount_genes += 1
					# low-memory option: don't store probes in memory: just print them
#					print '%s\t%s\t%s\t%i\t%i\t%s' %(probe,seq,strand,start,end,gene.id)
					# high memory option--allows further analysis and processing
					gene.probes.append( Probe(name,seq,strand,start,end) )
		sys.stderr.write('%i probes read, %i assigned in %i genes\n' %(pcount,pcount_assigned,pcount_genes) )

	def threePrime(self,nselect):
		# filter set of probes for each gene down to those nselect closest to 3' end of the gene
		for seq in self.genes:
#			print seq
			for strand in self.genes[seq]:
#				print strand
				for gene in self.genes[seq][strand]:
					probes = gene.probes
					probes.sort()
					nprobes=len(probes)
					if strand == '+':
						probes.reverse()
#						print seq,probes[0:nselect]
						gene.setProbes( probes[0:min(nselect,nprobes)] )
					elif strand == '-':
						gene.setProbes( probes[0:min(nselect,nprobes)] )

	def __str__(self):
		out = []
		for seq in self.genes:
#			print seq
			for strand in self.genes[seq]:
#				print strand
				for gene in self.genes[seq][strand]:
					out.append( str(gene) )
		return string.join(out)

	def printProbes(self):
		out = []
		for seq in self.genes:
			for strand in self.genes[seq]:
				for gene in self.genes[seq][strand]:
					for probe in gene.probes:
						out.append( '%s\t%s' %(str(probe),gene.id) )
		return string.join(out,'\n')

if __name__ == "__main__":
	p=OptionParser()
	p.add_option('-p','--probes')
	p.add_option('-g','--genes')
	p.add_option('-a','--antisense',action='store_true')
	p.add_option('--threeprime',type='int',default='3')
	opt,args=p.parse_args()
	app = App()
	app.run(opt,args)
	print app.printProbes()
