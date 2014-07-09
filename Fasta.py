#!/usr/bin/env python

__author__ = "Justin Ashworth"
# jashworth@systemsbiology.org
# Institute for Systems Biology

import os,re,string,sys
from optparse import OptionParser
from AshworthUtil import translate, rvs_comp_str
from InfoContent import *

class Fasta:
	sep='\t'
	nt = 'ACGTUacgtuNn'
	prot = 'ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy\.'
	re_nt = re.compile('^[%s]+$' %nt)
	re_notnt = re.compile('[^%s]' %nt)
	re_protein = re.compile('[%s]+' %prot)
	re_allprotein = re.compile('^[%s]+$' %prot)
	re_notprotein = re.compile('[^%s]' %prot)
	nttypes = ['dna','DNA','rna','RNA','nt','nucleotide']
	prottypes = ['protein','PROTEIN']
	testlength = 1000

	def __init__(self,seq='',name='',type=''):
		self.seq=seq
		self.name=name
		self.type=type
		self.clean()
		if type == '': self.guesstype()
		self.deepclean()

	def shortname(self,sep='|',field=-1):
		# names usually '|' separated with short/common name last (e.g. BLAST)
		fields=self.name.split(sep)
		name=re.split('[()\[\]]',fields[field])[0]
		return name.strip()

	def __getitem__(self,i):
		return self.seq[i]

	def __len__(self):
		return len(self.seq)

	def __eq__(self,other):
		# of self.name != other.name: return False
		return self.seq == other.seq

	def __cmp__(self,other):
		if self.name < other.name: return -1
		if self.name > other.name: return 1
		return 0

	def __hash__(self):
		return hash( self.name + self.seq )

	def subseq(self,start,end):
		l = len(self.seq)
		if start < 1 or start >= l or end < 1 or end > l:
			sys.stderr.write('warning: contig boundaries (%i,%i) invalid for %s (1,%i) [fixing]\n' %(start,end,self.name,l))
			if start < 1: start = 1
			if end < 1: end = 1
			if start >= l: start = l
			if end > l: end = l
		return Fasta( self.seq[start-1:end], '%s_%i-%i' %(self.name,start,end) )

	def rvscomp(self):
		self.seq = rvs_comp_str(self.seq)

	def guesstype(self):
		if self.seq == '': return
		l = min(self.testlength,len(self.seq))
		if self.re_nt.match(self.seq[:l]):
			self.type = self.nttypes[0]
		elif self.re_allprotein.match(self.seq[:l]):
			self.type = self.prottypes[0]
		else:
			sys.stderr.write('failed to guess type for %s\n' %self.name)

	def clean(self):
		self.seq = re.sub('[0-9 \t]','',self.seq)

	def deepclean(self):
		if self.type in self.nttypes:
			self.seq=self.re_notnt.sub('',self.seq)
		elif self.type in self.prottypes:
			self.seq=self.re_notprotein.sub('',self.seq)

	def GC(self):
		if not self.type in self.nttypes: return 0.0
		re_GC=re.compile('[GgCc]')
		return len( re_GC.findall(self.seq) )

	def GCc(self):
		if not self.type in self.nttypes: return 0.0
		return float(self.GC())/len(self.seq)

	def translate(self,frame=0):
		if frame < 0 or frame > 2:
			sys.stderr.write('ERROR: invalid frame selection %s: aborting translation (valid: 0,1,2)\n' %frame)
			return
		if not self.type in self.nttypes: return
		#self.seq = translate_simple(self.seq)
		translations = translate(self.seq)
		self.seq = string.join(translations[frame], '')
		self.type = 'protein'

	def summarize(self):
		summ = '%s:%s%s...%s%i' %(self.name,self.sep,self.seq[:10],self.sep,len(self.seq))
		if self.type in self.prottypes:
			summ += ' residues'
		elif self.type in self.nttypes:
			summ += ' bp'
			summ += '%s%0.2f GC' %(self.sep,self.GCc())
		return summ

	def fasta(self):
		return '>%s\n%s' %(self.name,self.seq)

	def csv(self):
		return '%s,%s' %(self.name,self.seq)

	def table(self,sep=None):
		if sep == None: sep = self.sep
		return '%s%s%s' %(self.name,sep,self.seq)

class SeqWindow:
	def __init__(self,name='',start=0,end=0,antisense=False):
		self.name=name
		self.start=start
		self.end=end
		self.antisense=antisense

class SeqWindows(dict):
	# file format expected:
	#five fields, tab/space-delimted: name source_seq start end strand
	#THAPSDRAFT_111 chr_1 400 799 +
	fwd = ['+','fwd','D']
	rvs = ['-','rvs','R']
	def __init__(self):
		dict.__init__(self)
		self.antisense=False
	def load(self,file):
		sys.stderr.write('reading subsequence definitions from %s\n' %file)
		for l in open(file):
			if l.startswith('#'): continue
			flds=l.strip().split()
			if len(flds) < 5: continue
			name = flds[0]
			source = flds[1]
			start = flds[2]
			end = flds[3]
			strand = flds[4]
			antisense = False
			if   strand in self.fwd: antisense = False
			elif strand in self.rvs: antisense = True
			else: sys.stderr.write('UNKNOWN DIRECTION ASSUMING fwd\n')
			if not self.has_key(source): self[source]=[]
			self[source].append( SeqWindow( name, int(start), int(end), antisense ))

class FastaSeqs:
	op=OptionParser()
	op.add_option('-i','--infocontent',action='store_true')
	op.add_option('--filter',help='regex filter for sequence names')
	op.add_option('--ids',help='file containing list (single column) of specific ids to consider')
	op.add_option('-O','--overwrite',action='store_true')
	op.add_option('-r','--raw',help='just output space-delimited raw sequences',action='store_true')
	op.add_option('--rawfiles',help='write each fasta sequence to a separate, unformatted file',action='store_true')
	op.add_option('-s','--summarize',action='store_true')
	op.add_option('-S','--separate',action='store_true')
	op.add_option('--subseqs')
	op.add_option('--upstream')
	op.add_option('-t','--translate',action='store_true')
	op.add_option('-f','--frame',type='int',default=0,help='reading frame offset for translation')
	op.add_option('-u','--unique',action='store_true')
	op.add_option('--rvscomp',action='store_true')
	op.add_option('--csv',action='store_true')

	def __init__(self):
		self.seqs={}
		# IMPORTANT: LOOP AND OUTPUT FUNCTIONS USE THIS INDEX (so be careful when changing 'seqs')
		self.order=[]

	def go(self,opt,args):
		self.loadseqs(args)
		self.process(opt)
		if len(self.seqs)==0: return "no sequences"
		elif opt.raw: return self.rawseqs()
		elif opt.rawfiles: return self.rawfiles()
		elif opt.separate: return self.separate(opt.overwrite)
		elif opt.summarize: return self.summarize()
		elif opt.infocontent: return self.infocontent()
		elif opt.csv: return self.csv()
		else: return str(self)

	def loadseqs(self,files):
		for file in files:
			if not file or not os.path.exists(file):
				print 'skipping invalid sequence file %s' %file
				continue
			seqname=''
			seq=[]
			sys.stderr.write('reading sequence from %s\n' %file)
			for line in open(file):
				if line.startswith('>'):
					if seqname != '' and len(seq)>0:
						if self.seqs.has_key(seqname):
							sys.stderr.write('warning: replacing exisiting seq with same name %s\n' %seqname)
						self.seqs[seqname] = Fasta( string.join(seq,''), seqname )
						self.order.append(seqname)
						seq=[]
					seqname=re.sub('>','',line.strip()).strip()
				else:
					seq.append( line.strip() )
			# last hanging seq
			if seqname != '' and len(seq)>0:
				if self.seqs.has_key(seqname):
					sys.stderr.write('warning: replacing exisiting seq with same name %s\n' %seqname)
				self.seqs[seqname] = Fasta( string.join(seq,''), seqname )
				self.order.append(seqname)

	def process(self,opt):
		if opt.subseqs:
			windows=SeqWindows()
			windows.load(opt.subseqs)
			subseqs = self.subseqs(windows)
			# deep copy/replace
			self.seqs = subseqs.seqs
			self.order = subseqs.order
		if opt.translate: self.translate(opt.frame)
		if opt.filter: self.filter(opt)
		if opt.ids: self.these_ids(opt.ids)
		if opt.unique: self.unique()
		if opt.rvscomp: self.rvscomp()

	def translate(self,frame=0):
		for seq in self.seqs.values(): seq.translate(frame)

	def filter(self,opt):
		filteredseqs = {}
		neworder = []
		re_filter = re.compile(opt.filter)
		for name,seq in self.seqs.items():
			if re_filter.search(name):
				sys.stderr.write('%s\n' %name)
				filteredseqs[name] = seq
				neworder.append(name)
		if len(filteredseqs) == 0:
			sys.stderr.write('oops, no seq names passed filter %s!\n' %opt.filter)
		self.seqs = filteredseqs
		self.order = neworder

	def these_ids(self,ids_file):
		if not os.path.exists(ids_file):
			sys.stderr.write('ids file %s not found' %ids_file)
			return
		ids = [ l.strip() for l in open(ids_file) ]
		keep_seqs = {}
		neworder = []
		for name,seq in self.seqs.items():
			if not name in ids: continue
			keep_seqs[name]=seq
			neworder.append(name)
		self.seqs = keep_seqs
		self.order = neworder

	def unique(self):
		uniqueseqs = {}
		neworder = []
		for seq in self.seqs.values():
			unique = True
			for useq in uniqueseqs.values():
				if seq == useq:
					unique = False
					break
			uniqueseqs[seq.name] = seq
			neworder.append(seq.name)
		self.seqs = uniqueseqs

	def rvscomp(self):
		for seq in self.seqs.values():
			seq.rvscomp()

	def rawseqs(self):
		return string.join( [self.seqs[key].seq for key in self.order],' ')

	def subseqs(self,windows):
		subseqs = FastaSeqs()
		for seq,windows in windows.items():
			if not self.seqs.has_key(seq): print 'no seq %s loaded' %seq; continue
			for window in windows:
				subseq = self.seqs[seq].subseq(window.start,window.end)
				subseq.name = window.name
				if window.antisense: subseq.rvscomp()
				subseqs.seqs[subseq.name] = subseq
				subseqs.order.append(subseq.name)
		return subseqs

	def separate(self,overwrite=False):
		for key in self.order:
			seq = self.seqs[key]
			fname='%s.fa'%seq.shortname()
			if os.path.exists(fname):
				if overwrite: print 'overwriting file %s' %fname
				else:
					print '%s exists, skipping' %fname
					continue
			f=open(fname,'w')
			f.write('%s\n' %seq.fasta())
			f.close()
		return self.summarize()

	def rawfiles(self):
		for key in self.order:
			seq = self.seqs[key]
			f=open('%s.raw'%seq.shortname(),'w')
			f.write('%s\n'%seq.seq.lower())
			f.close()
		return self.summarize()

	def csv(self):
		out = []
		for key in self.order:
			seq = self.seqs[key]
			out.append( seq.csv() )
		return string.join( out, '\n' )

	def __str__(self):
		out = []
		for key in self.order:
			seq = self.seqs[key]
			out.append( seq.fasta() )
		return string.join( out, '\n' )

	# Schneider/Stormo/(Shannon) information content
	def infocontent(self):
		# Information content is formally the Shannon sequence informational entropy of the query sequence vs. the background sequence. However, the background entropy is usually maximal, and equal to the sum(P*logP) where all P are equal. This value is 2 for log base 2, or 1 for log base 4.
		# IC is calculated in base 2 in range [0,2] for historical reasons.
#		ic = InfoContent(2,2)
		# However, for DNA sequences it can also be calculated in [0,1] with base 4 (where 1 = maximal background entropy for DNA sequence using base 4)
		ic = InfoContent(4,1)
		rawseqs = []
		for seq in self.seqs.values():
			rawseqs.append(seq.seq)
		values = ic.run(rawseqs)
		return string.join(['%.3f'%v for v in values],'\n')

	def summarize(self):
		return string.join([ self.seqs[key].summarize() for key in self.order], '\n' )

	def named_list(self):
		seqs = []
		for key in self.order:
			seq = self.seqs[key]
			seqs.append( (key,seq) )
		return seqs

	def simple_list(self):
		seqs = []
		for s in self.seqs.values():
			seqs.append(s.seq)
		return seqs

################################################################################

if __name__ == "__main__":
	op=FastaSeqs.op
	opt,args=op.parse_args()

	app=FastaSeqs()
	print app.go(opt,args)
