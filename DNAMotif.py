#!/usr/bin/env python
import os,re,string,sys,math
from optparse import OptionParser
from Fasta import *
import random

'''
DESCRIPTION
Input: Fasta file of DNA sequences of equal length
Output: DNA 'motif' count, probability, or log-odds matrices (such as for TRANSFAC or MEME/MAST)
Can be used instead of transfac2meme
To do: merge/interface with InfoContent class
'''

__author__ = "Justin Ashworth"

RandPos = { 'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25 }

class DNAMotif:
	bases = ['A','C','G','T']
	bg_default = { 'A':0.25, 'C':0.25, 'G':0.25, 'T':0.25 }
	# zero-substitution allows log-odds calculations, has similar effect to using pseudocounts
	ZERO_PROB = 1e-2
	sep = '\t'
	#sep = ' '
	def __init__(self,width=0,name='',pseudocounts=0,data_type='',bg=None):
		self.set_width(width)
		self.name = name
		self.pseudocounts = pseudocounts
		self.data_type = data_type
		self.bg = bg
		if not self.bg: self.bg = self.bg_default
		self.score = None

	def set_width(self,width):
		self.width = width
		self.init_matrix()

	def init_matrix(self):
		self.matrix = {
			'counts' : [],
			'probs' : [],
			'logodds' : [],
		}
		for i in range(self.width):
			self.matrix['counts' ].append({ 'A':0,   'C':0,   'G':0,   'T':0 })
			self.matrix['probs'  ].append({ 'A':0.0, 'C':0.0, 'G':0.0, 'T':0.0 })
			self.matrix['logodds'].append({ 'A':0.0, 'C':0.0, 'G':0.0, 'T':0.0 })

	def truncate(self,width):
		width = min(self.width,width)
		for type in self.matrix:
			self.matrix[type] = self.matrix[type][:width]
		self.width = width

	def weaken(self,value):
		#for in range(len(self.matrix['probs'])):
		#	for b,val in self.matrix['probs'][i].items():
		#		self.matrix['probs'][i][b] += value
		for pos in self.matrix['probs']:
			for b,val in pos.items():
				pos[b] += value
		self.normalize_probs()

	def pad_left(self,pad):
		w = len(self.matrix['probs'])
		newmat = []
		for i in range(pad):
#			print i
			newmat.append(RandPos)
		for i in range(pad,w+pad):
#			print i
			newmat.append(self.matrix['probs'][i-pad])
		self.set_width(w+pad)
		self.matrix['probs'] = newmat
#		print 'padded left by %i: %s' %(pad, self.consensus())

	def pad_right(self,pad):
		for i in range(pad):
			self.matrix['probs'].append(RandPos)
		self.width=len(self.matrix['probs'])
#		print 'padded right by %i: %s' %(pad, self.consensus())

	def consensus(self):
		cons = []
		for pos in self.matrix['probs']:
			rand = True
			topb = 'X'
			topv = 0
			for b,p in pos.items():
				if not p==0.25: rand=False
				if p>topv:
					topv=p
					topb=b
			if rand: cons.append('n')
			elif topb=='X': cons.append('X')
			else: cons.append(topb)
		return string.join(cons,'')

	def matrix_from_lines(self,lines,type):
		self.set_width(len(lines))
		self.data_type = type
		for i in range(len(lines)):
			l = lines[i]
			vals = l.strip().split()
			#print vals
			for b in range(len(self.bases)):
				self.matrix[type][i][self.bases[b]]=float(vals[b])
		self.normalize_probs()
		print 'read %s' %self.consensus()

	def load_motif_file(self,fname):
		if not os.path.exists(fname): raise IOError('file %s not found' %fname)
		format=''
		for line in open(fname):
			if line.startswith('#name='): self.name = line.strip().split('=')[-1]
			if line.startswith('#format='): format = line.strip().split('=')[-1]
			if line.startswith('#type='): self.data_type = line.strip().split('=')[-1]
		if format == 'simple':
			self.load_motif_simple(fname,self.data_type)
		else:
			sys.stderr.write('loading motif from file with unspecified format')
			self.load_motif_simple(fname)

	def load_motif_simple(self,fname,data_type):
		sys.stderr.write('loading simple matrix from: %s as type: %s\n' %(fname,data_type) )
		colnames = []
		temp_matrix = []
		colnames = ['a','c','g','t']
		for line in open(fname):
			if line.startswith('#'): continue
			if len(line) == 0: continue
			if line.startswith('key'):
				colnames = line.strip().split()[1:]
				continue
			fields = line.strip().split()
			if not len(fields) == len(colnames):
				raise RuntimeError('number of fields doesn\'t match number of colnames')
			rowdict = {}
			for i in range(len(fields)):
				rowdict[colnames[i]] = float(fields[i])
			temp_matrix.append(rowdict)
		if data_type in self.matrix:
			self.set_width(len(temp_matrix))
			self.matrix[data_type] = temp_matrix
		else: raise RuntimeError('unknown/unspecified matrix type %s' %data_type)

	def read_bgfile(self,file):
		bg = { 'A':-1, 'C':-1, 'G':-1, 'T':-1 }
		for line in open(file):
			words = line.strip().split()
			if line.startswith('#'): continue
			if len(words) != 2: raise IOError('bad bgfile')
			base = words[0]
			if not base in self.bases: raise IOError('bad bgfile')
			frac = float(words[1])
			if frac < 0 or frac > 1: raise RuntimeError('bad bgfile')
			bg[base] = frac
			total = 0
			for base,frac in bg.items():
				if frac == -1: RuntimeError('bad bgfile')
				total += frac
			if total != 1.0: RuntimeError('bad bgfile')
		self.bg = bg
		return(bg)

	def make_counts_matrix(self,sequences=[]):
		self.set_width( len(sequences[0]) )
		for s in sequences:
			#print s
			l = len(s)
			if l != self.width: raise RuntimeError('length mismatch')
			for i in range(l):
				base = s[i].upper()
				if not base in self.bases: raise RuntimeError('wrong alphabet')
				self.matrix['counts'][i][base] += 1

	def probs_from_seqs(self,seqs):
		self.make_counts_matrix(seqs)
		self.make_probs_matrix()

	def make_probs_matrix(self):
		# (meaningless without filled counts matrix)
		for i in range(self.width):
			total = 0
			for base in self.bases:
				total += self.matrix['counts'][i][base]
			for base in self.bases:
				# as in transfac2meme
				self.matrix['probs'][i][base] = float( self.matrix['counts'][i][base] + self.pseudocounts * self.bg[base] ) / ( total + self.pseudocounts )

	def normalize_probs(self):
		for pos in self.matrix['probs']:
			tot = sum(pos.values())
			for b,prob in pos.items():
				pos[b] = float(prob)/tot

	def make_logodds_matrix(self):
		# (meaningless without filled probs matrix)
		if self.data_type == 'counts': self.make_probs_matrix()
		for i in range(self.width):
			for base in self.bases:
				# as in transfac2meme
				#log( $motif{$i_base, $i_motif} / $bg{$bases[$i_base]} ) /log(2.0)
				prob = self.matrix['probs'][i][base]
				if prob == 0.0: prob = self.ZERO_PROB
				self.matrix['logodds'][i][base] = math.log( float(prob) / self.bg[base] ) / math.log(2.0)

	def __str__(self):
		return self.output_default('probs')

	def __len__(self):
		return len(self.matrix['probs'])

	def same_probs(self,other):
		p1 = self.matrix['probs']
		p2 = other.matrix['probs']
		for i in range(len(p1)):
			for b in self.bases:
				v1 = float(p1[i][b])
				v2 = float(p2[i][b])
				if abs(v2-v1) > 1e-10:
#					print '%i,%s,diff %g' %(i,b,v2-v1)
					return False
		return True

	def output(self,type='',name=None):
		if not name: name = self.name
		if type == 'MEME': return self.output_MEME(name)
		elif type == 'TRANSFAC': return self.output_TRANSFAC(name)
		else: return self.output_default(type,name)

	def output_default(self,type='probs',name=None):
		out = []
		if not self.data_type in self.matrix: raise RuntimeError('no data_type set')
		if type == 'probs' and self.data_type=='counts': self.make_probs_matrix()
		out.append('#name=%s' %name)
		out.append('#format=simple')
		out.append('#type=%s'%type)
		out.append('key%s%s' %(self.sep,string.join(self.bases,self.sep)))
		for i in range(self.width):
			line = str(i)
			for base in self.bases:
				value = ''
				if type=='counts' and self.data_type=='counts': value = '%i' %self.matrix['counts'][i][base]
				elif type=='probs': value = '%.4f' %self.matrix['probs'][i][base]
				else: raise RuntimeError('undefined output request')
				line += '%s%s' %(self.sep,value)
			out.append(line)
		return string.join(out,'\n') + '\n'

	def output_TRANSFAC(self,name):
		# transfac2meme chokes on this as currently formatted(?)
		out = []
		out.append('ID%s%s' %(self.sep,name))
		out.append('PO%s%s' %(self.sep,string.join(self.bases)))
		for i in range(self.width):
			line = str(i+1)
			for base in self.bases:
				line += '%s%i' %(self.sep,self.matrix['counts'][i][base])
			line += '%s-' %self.sep
			out.append(line)
		out.append('//')
		return string.join(out,'\n') + '\n'

	def output_MEME_header(self):
		out = []
		out.append( 'MEME version 4\n' )
		out.append( 'ALPHABET= %s\n' %(string.join(self.bases,'')) )
		out.append( 'strands: + -\n' )
		out.append( 'Background letter frequencies (from' )
		out.append( string.join( [ '%s %f' %(base,self.bg[base]) for base in self.bases ] ) )
		return string.join(out,'\n') + '\n'

	def output_MEME(self,name=None,logodds=False):
#		name=re.sub('\.','_',name)
#		name=re.sub('\.','',name)
		if not name: name = self.name
		out = []
		out.append( 'MOTIF %s' %name )
		if logodds: out.append( self.output_MEME_logodds() )
#		out.append( 'MOTIF %s' %name )
		out.append( self.output_MEME_probs() )
		return string.join(out,'\n') + '\n'

	def output_MEME_logodds(self):
		# background-corrected, psuedocount log-odds scores for MEME/MAST
		out = []
		self.make_logodds_matrix()
		alength = len(self.bases)
		out.append( 'log-odds matrix: alength= %i w= %i ' %(alength,self.width) )
		for i in range(self.width):
			line = []
			for base in self.bases:
				line.append('%f' %self.matrix['logodds'][i][base])
			out.append(string.join(line,self.sep))
		return string.join(out,'\n') + '\n'

	def output_MEME_probs(self):
		out = []
		if self.matrix['probs']==[]: self.make_probs_matrix()
		alength = len(self.bases)
		out.append( 'letter-probability matrix: alength= %i w= %i nsites= 1 E= 0 ' %(alength,self.width) )
		for i in range(self.width):
			line = []
			for base in self.bases:
				line.append('%f' %self.matrix['probs'][i][base])
			out.append(string.join(line,self.sep))
		return string.join(out,'\n') + '\n'

def mix_motifs_average(mot1,mot2,type='probs'):
	mat1 = mot1.matrix[type]
	mat2 = mot2.matrix[type]
	if not len(mat1) == len(mat2):
		sys.stderr.write('ERROR: can\'t mix motifs of different lengths\n')
		return None
	child = DNAMotif(width=len(mat1))
	for i in range(len(mat1)):
		for b in mot1.bases:
			child.matrix[type][i][b] = float(mat1[i][b] + mat2[i][b]) / 2
	child.normalize_probs()
	child.data_type = 'probs'
	return child

def mix_motifs_recombine(mot1,mot2,type='probs'):
	mat1 = mot1.matrix[type]
	mat2 = mot2.matrix[type]
	l1 = len(mat1)
	l2 = len(mat2)
	l = min(l1,l2)
	child = DNAMotif(width=l)
	for i in range(l):
		# take position at random from a source mot, biased by score
		s1 = mot1.score
		s2 = mot2.score
		rnd = random.uniform(0,s1+s2)
		if rnd < s1:
			for b in mot1.bases: child.matrix[type][i][b] = mat1[i][b]
		else:
			for b in mot1.bases: child.matrix[type][i][b] = mat2[i][b]
	child.normalize_probs()
	child.data_type = 'probs'
	return child

def mix_motifs(mot1,mot2,method='average',type='probs'):
	if method=='average':
		return mix_motifs_average(mot1,mot2,type)
	if method=='recombine':
		return mix_motifs_recombine(mot1,mot2,type)

def center_gapped_motif(mot,gap):
	newmot = DNAMotif()
	w = len(mot.matrix['probs'])
	mid = int(w/2)
	gapi = range(mid,mid+gap)
	newmot.set_width(w+gap)
#	print w, mid
#	print gapi
	for i in range(mid):
#		print i
		newmot.matrix['probs'][i] = mot.matrix['probs'][i]
	for i in gapi:
#		print 'gap',i
		newmot.matrix['probs'][i] = RandPos
	for i in range(mid+gap,w+gap):
#		print i
		newmot.matrix['probs'][i] = mot.matrix['probs'][i-gap]
	print 'created center-gapped motif: %s' %newmot.consensus()
	newmot.name = '%s_cg%i' %(mot.name,gap)
	newmot.data_type = 'probs'
	return newmot

class DNAMotifs:
	def __init__(self):
		self.motifs = []

	def run(self,opt,args):
#		if len(args) < 1: raise IOError('no arguments')
		self.opt = opt
		if opt.motif_file:
			motif = DNAMotif()
			if opt.bg:
				if not os.path.exists(opt.bg): RuntimeError('bg file not found')
				motif.read_bgfile(opt.bg)
			motif.load_motif_file(opt.motif_file)
			self.motifs.append(motif)
		else: self.load_fasta_files(args,opt)

	def load_fasta_files(self,files,opt=None):
		if opt: self.opt = opt
		bg=DNAMotif.bg_default
		if opt.bg:
			if not os.path.exists(opt.bg): RuntimeError('bg file not found')
			bg=DNAMotif().read_bgfile(opt.bg)

		for fname in files:
			if not os.path.exists(fname): raise IOError('file %s not found' %fname)
			dnaseqs = FastaSeqs()
			dnaseqs.loadseqs([fname])
			seqlen = 0
			seqs = dnaseqs.seqs.values()

			if opt.start or opt.end and opt.end != 0:
				if opt.start <= 0 or opt.start >= seqlen: RuntimeError('bad start')
				if opt.end <= 0 or opt.end >= seqlen: RuntimeError('bad end')
				seqs = [ seq.subseq(opt.start,opt.end) for seq in seqs ]
#				print string.join([s.seq for s in seqs],'\n')

			for seq in seqs:
				l = len(seq)
				if seqlen == 0: seqlen = l
				elif l != seqlen: raise RuntimeError('length mismatch')

			motif = DNAMotif(seqlen,fname,pseudocounts=opt.pseudocounts,data_type='counts',bg=bg)

			motif.make_counts_matrix( [seq.seq for seq in seqs] )
			self.motifs.append(motif)

	def __str__(self):
		out = []
		for m in self.motifs:
			out.append(str(m))
		return string.join(out,'\n') + '\n'

	def output(self,type='counts'):
		out = []
		if type == 'MEME': out.append(self.motifs[0].output_MEME_header())
		for m in self.motifs:
			name = m.name
			if opt.outprefix: name = opt.outprefix
			out.append(m.output(type=type,name=name))
		out = string.join(out,'\n') + '\n'
		if self.opt.outprefix:
			fname = '%s.%s' %(opt.outprefix,type)
			f=open(fname,'w')
			f.write(out)
			f.close()
			return 'output %s motif to %s' %(type,fname)
		else: return out

def readMEME_from_lines(lines,Eval_cutoff=1e100):
	mots = []
	record = False
	matlines = []
	re_eval = re.compile('E= ([0-9+-.e]+)')
	re_name = re.compile('MOTIF\s+([^\s]+)')
	Eval = 1e100
	for line in lines:
		namematch = re_name.search(line)
		if namematch:
			name = namematch.groups()[0]
			print 'reading %s' %name
		if line.startswith('letter-probability matrix'):
			Eval = re_eval.search(line)
			if Eval: Eval = float( Eval.groups()[0] )
			else: Eval = None
			record = True
			continue
		if line.startswith('------') or len(line.strip().split())<4:
			if record:
				if Eval <= Eval_cutoff or not Eval:
					mot = DNAMotif(name='%s_%i'%(name,len(mots)))
					mot.matrix_from_lines(matlines,'probs')
					mots.append(mot)
				# reset
				Eval = 1e100
				matlines = []
			record = False
		if record: matlines.append(line)
	if record:
		if Eval <= Eval_cutoff or not Eval:
			mot = DNAMotif(name='%s_%i'%(name,len(mots)))
			mot.matrix_from_lines(matlines,'probs')
			mots.append(mot)
	return(mots)

def readMEME(path,Eval_cutoff=1e100):
	return readMEME_from_lines(open(path).readlines())

if __name__ == "__main__":
	p=OptionParser()
	p.add_option('--bg',help='background base frequencies')
	p.add_option('--end',type='int',default=0,help='motif endpoint in sequence')
	p.add_option('--MEME',action='store_true')
	p.add_option('-m','--motif_file')
	p.add_option('-p','--probs',action='store_true')
	p.add_option('-o','--outprefix')
	p.add_option('--pseudocounts',type='int',default='0')
	p.add_option('--start',type='int',default=0,help='motif starting point in sequence')
	p.add_option('-t','--TRANSFAC',action='store_true')
	opt,args=p.parse_args()
	app = DNAMotifs()
	app.run(opt,args)
	if opt.MEME: print app.output('MEME')
	elif opt.probs: print app.output('probs')
	elif opt.TRANSFAC: print app.output('TRANSFAC')
	else: print app.output()
