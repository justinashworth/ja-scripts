#!/usr/bin/env python
__author__ = "Justin Ashworth"

import commands, re, sys, string, os
from math import sqrt

# old python versions do not have the 'sorted' built-in
try:
	if sorted: pass
except:
	def sorted( list, func = None ):
		list.sort( func )
		return list

try:
	if reversed: pass
except:
	def reversed(old):
		new = []
		l = len(old)
		for i in range(l):
			new.append(old[l-i-1])
		return new

# older versions of python may not have the 'set' function
def simple_set_fxn(list):
	set = []
	for thing in list:
		if not thing in set: set.append( thing )
	return set
try:
	if set: pass
except:
	def set(list):
		return simple_set_fxn(list)

def jglob( regex, positive = True, dir = '.' ):
	match = re.compile(regex)
	prefix = ''
	if dir != '.': prefix = dir+'/'
	if positive:
		return sorted( [ prefix+file for file in os.listdir(dir) \
		                 if not os.path.isdir(file) \
		                 and not file.startswith('.') \
		                 and match.search(file) ] )
	return sorted( [ prefix+file for file in os.listdir(dir) \
	                 if not os.path.isdir(file) \
	                 and not file.startswith('.') \
	                 and not match.search(file) ] )

def recurse_dirs(dir,func):
	for nextdir in os.listdir(dir):
		if nextdir == '.': continue
		if nextdir == '..': continue
		if os.path.isdir(nextdir):
			recurse_dirs(nextdir,func)
	func(dir)

re_gz = re.compile('gz')
re_letters = re.compile('[a-zA-Z]')
re_numbers = re.compile('[0-9]')

aa1 = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
if not len(aa1) == 20:
	print 'error!'; sys.exit

aa3 = ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR']
if not len(aa3) == 20:
	print 'error!'; sys.exit
amino_acids = aa3

oneletter = {'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G',
             'HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N',
             'PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V',
             'TRP':'W','TYR':'Y',
						 'ADE':'a','CYT':'c','GUA':'g','THY':'t',
						 '  A':'a','  C':'c','  G':'g','  T':'t',
						 'A':'a','C':'c','G':'g','T':'t',
						 'DA':'a','DC':'c','DG':'g','DT':'t', # new pdb spec
						 ' DA':'a',' DC':'c',' DG':'g',' DT':'t'} # new pdb spec

def get_oneletter(aa3):
	if oneletter.has_key(aa3): return oneletter[aa3]
	return aa3.strip()[0]

threeletter = {}
threeletter = {'X':'STP'}
for three,one in oneletter.items(): threeletter[one] = three
threeletter['a'] = 'ADE'
threeletter['c'] = 'CYT'
threeletter['g'] = 'GUA'
threeletter['t'] = 'THY'

def get_threeletter(letter):
	if threeletter.has_key(letter): return threeletter[letter]
	return 'UNK'

pdb_residues = amino_acids + [' DA',' DC',' DG',' DT']

bases = ['a','c','g','t']
BASES = ['A','C','G','T']

def isDNA(query):
	if query in bases: return True
	if query in BASES: return True
	if query in dnares: return True
	return False

comp = {
	'A':'T',
	'C':'G',
	'G':'C',
	'T':'A',
	'a':'t',
	'c':'g',
	'g':'c',
	't':'a',
#	'n':'n', # use compbase if other characters are expected
#	'N':'N',
}

# safer
def compbase(char):
	try: return comp[char]
	except: return char

# redundant with oneletter?
dnares = {
	'  A':'A', 'ADE':'A', ' DA':'A',
	'  C':'C', 'CYT':'C', ' DC':'C',
	'  G':'G', 'GUA':'G', ' DG':'G',
	'  T':'T', 'THY':'T', ' DT':'T',
}
nucleic_acids = dnares.keys()

codons = {
	'A' : ['GCA','GCC','GCG','GCT'],
	'C' : ['TGC','TGT'],
	'D' : ['GAC','GAT'],
	'E' : ['GAA','GAG'],
	'F' : ['TTC','TTT'],
	'G' : ['GGA','GGC','GGG','GGT'],
	'H' : ['CAC','CAT'],
	'I' : ['ATA','ATC','ATT'],
	'K' : ['AAA','AAG'],
	'L' : ['CTA','CTC','CTG','CTT','TTA','TTG'],
	'M' : ['ATG'],
	'N' : ['AAC','AAT'],
	'P' : ['CCA','CCC','CCG','CCT'],
	'Q' : ['CAA','CAG',],
	'R' : ['CGA','CGC','CGG','CGT','AGA','AGG'],
	'S' : ['AGC','AGT','TCA','TCC','TCG','TCT'],
	'T' : ['ACA','ACC','ACG','ACT'],
	'V' : ['GTA','GTC','GTG','GTT'],
	'W' : ['TGG'],
	'Y' : ['TAC','TAT'],
	'.' : ['TAA','TAG','TGA'],
}

degen_nucs = {
	'R' : set(['A','G']),
	'Y' : set(['C','T']),
	'M' : set(['A','C']),
	'K' : set(['G','T']),
	'S' : set(['C','G']),
	'W' : set(['A','T']),
	'B' : set(['C','G','T']),
	'D' : set(['A','G','T']),
	'H' : set(['A','C','T']),
	'V' : set(['A','C','G']),
	'N' : set(['A','C','G','T']),
}

def degen_code_for_nucs(nucs):
	nucset = set( [ n for n in nucs ] )
	for code,codenucs in degen_nucs.items():
		if nucset == codenucs: return code
	return '?'

def pdb_to_seq(filename,template=None):
	chains = {}
	re_calpha = re.compile('^ATOM\s+[0-9]+  CA [ A]([A-Z]+) ([A-Z])([0-9 -]+) ')
	re_c1star = re.compile('^ATOM\s+[0-9]+  C1[\*\']\s+([A-Z]+) ([A-Z])([0-9 -]+) ')

	pdb = None
	if re.search('gz$',filename):
		pdb = os.popen('gunzip -c %s' %filename)
	else: pdb = file(filename,'r')

	if template and os.path.exists(template):
		print 'template file %s' %template
		chains = pdb_to_seq(template) # semi-recursive

	for line in pdb:
		match = re_calpha.search(line) # look for protein c-alpha
		if match != None:
			resn,chain,resi = match.groups()
			if not chains.has_key('prot'): chains['prot'] = {}
			if not chains['prot'].has_key(chain): chains['prot'][chain] = Chain()
			letter = 'X'
			try: letter = oneletter[resn]
			except: pass
			chains['prot'][chain].add( ChainPos(chain,resi), Residue(letter) )
			continue
		match = re_c1star.search(line) # look for DNA c1* or c1'
		if match != None:
			resn,chain,resi = match.groups()
			if not chains.has_key('dna'): chains['dna'] = {}
			if not chains['dna'].has_key(chain): chains['dna'][chain] = Chain()
			letter = 'X'
			try: letter = oneletter[resn]
			except: pass
#			print line,chain,resi,letter
			chains['dna'][chain].add( ChainPos(chain,resi), Residue(letter) )
	return chains

def pdb_seq_dict(file):
	chains = {}
	re_calpha = re.compile('^ATOM\s+[0-9]+  CA [ A]([A-Z]+) ([A-Z])([ \-0-9]{4})')
	re_c1star = re.compile('^ATOM\s+[0-9]+  C1[\*\']\s+([A-Z]+) ([A-Z])([ \-0-9]{4})')
	for line in open(file):
		match = re_calpha.search(line) # look for protein c-alpha
		if match != None:
			(resn,chain,resi) = match.groups()
			if not chains.has_key(chain): chains[chain] = {}
			chains[chain][int(resi)] = oneletter[resn]
			continue
		match = re_c1star.search(line) # look for DNA c1*
		if match != None:
			(resn,chain,resi) = match.groups()
			if not chains.has_key(chain): chains[chain] = {}
			chains[chain][int(resi)] = oneletter[resn]
	return chains

# grab information related to protein-DNA design parameters out of the rosetta filename, expecting some ambiguity as to the order and positioning of keys in the filename
def rosetta_filename(filename):
	ids = {}
	key_re = [
		# order matters
		( 'queueindex', re.compile('^[0-9]+$') ),
		( 'pdb',        re.compile('^[a-zA-Z0-9]{4}$') ),
		( 'chain',      re.compile('^[A-Z]$') ),
		( 'pos',        re.compile('^[0-9-]+$') ),
		( 'seq',        re.compile('^[acgtACGTm]+$') ),
		( 'newseq',     re.compile('^[A-Z]\.[0-9-]+\.[a-zA-Z]+') ),
		( 'index',      re.compile('^[0-9]{4}$') )
	]

	pathsplit = filename.split('/')
	for pathfield in pathsplit:
#		undsplit = pathfield.split('.')[0].split('_')
		undsplit = pathfield.split('_')
		for undfield in undsplit:
			for key,regex in key_re:
				if not ids.has_key(key):
					if regex.match(undfield):
						ids[key] = undfield
						break

	ids['gz'] = False
	if filename.split('.')[-1] == 'gz': ids['gz'] = True

	# string->integer conversions
	for key in ['queueindex','pos','index']:
		if ids.has_key(key): ids[key] = int( ids[key] )

	prefix = []
#	for id in ['pdb','pos','seq','index']:
	for id in ['pdb','pos','seq']:
#	for id in ['pdb','pos','seq','queueindex']:
		if ids.has_key(id): prefix.append( str(ids[id]) )

	prefix = string.join( prefix, '_' )
#	print prefix, ids
	return prefix, ids

def rosetta_filenames_sort(a,b):
	keys = [ 'pdb', 'chain', 'pos', 'seq', ]
	ids_a = rosetta_filename(a)[1]
	ids_b = rosetta_filename(b)[1]
	for key in keys:
		if not ids_a.has_key(key) or not ids_b.has_key(key): continue
		if ids_a[key] < ids_b[key]: return -1
		if ids_a[key] > ids_b[key]: return 1
	return 0

def is_hydrogen(atomname):
# checks whether first letter is 'H' (all rosetta hydrogens)
	for l in atomname:
		if re_numbers.match(l) or l == ' ': continue
		if l == 'H': return True
		else: return False

bbatoms = ['N','CA','C','O','CB']
#bbatoms = ['N','CA','C','O']

def skip_atom(atomname):
	if is_hydrogen(atomname): return True
	if atomname in bbatoms: return True
	return False

translation_code = {}
n = 0
for aa,codonlist in codons.items():
	n += len(codonlist)
	for codon in codonlist:
		translation_code[codon] = aa

def get_codons(aa):
	if codons.has_key(aa): return codons[aa]
	return ['???']

def translate_codon(codon):
	if translation_code.has_key( codon.upper() ): return translation_code[ codon.upper() ]
	return "X"

def stringf_left( string, length ):
	slen = len(string)
	if slen > length: return string
	for i in range( length - slen ): string += ' '
	return string

def regex_remove( regex, string ):
	if regex.search( string ): return False
	return True

def regex_filter( regex, string ):
	return regex.sub( '', string )

# figure out the next dna sequence after 'seq', according to rosetta ordering

next_rosetta_base = { 'g':'a', 'a':'c', 'c':'t', 't':'g' }

def next_seq( seq ):
	i = len(seq)-1
	next = []
	for j in range(i+1): next.append('')
	change = 1
	while i >= 0:
		next[i] = seq[i]
		if change == 1: next[i] = next_rosetta_base[seq[i]]
		if seq[i] != 't': change = 0
		i -= 1
	return string.join(next,'')

nucs = ['a','c','g','t','n']

def lead_zero_string( length, number ):
	string = str(number)
	if len(string) > length:
		print 'number greater than length!'
		sys.exit()
	while len(string) < length:
		string = '0' + string
	return string

class Position:
	def __init__(self,chain=' ',index=0,type='N'):
		self.chain = chain
		self.index = index
		self.type = type
	def __str__(self):
		return '%s.%i.%s' %(self.chain,self.index,self.type)

# recursive sequence-generators
def all_combinations( pos, length, seq, seqs ):
	if seq == []:
		for i in range(length): seq.append( '' )
	for base in bases:
		seq[pos] = base
		if pos == length-1:
			seqs.append( string.join( seq, '' ) )
		else: all_combinations( pos+1, length, seq, seqs )

def all_combinations_gen( pos, length, seq, seqs, choices ):
	if seq == []:
		for i in range(length): seq.append( '' )
	for choice in choices[pos]:
		seq[pos] = choice
		if pos == length-1:
			seqs.append( string.join( seq, '' ) )
		else: all_combinations_gen( pos+1, length, seq, seqs, choices )

def all_combinations_without( pos, length, seq, seqs, nativeseq ):
	if seq == []:
		for i in range(length): seq.append( '' )
	for base in bases:
		if base == nativeseq[pos]: continue
		seq[pos] = base
		if pos == length-1:
			seqs.append( string.join( seq, '' ) )
		else: all_combinations_without( pos+1, length, seq, seqs, nativeseq )

def only_combos_with( base, ipos, seqs ):
#	def has_the_base( seq ):
#		return seq[ipos] == base
#	return filter( has_the_base, seqs )
	return [ seq for seq in seqs if seq[ipos] == base ]

def permute_positions( index, seedseq, seq, seqs ):
	# the elements of seqs are of type Position (above)
	# index: current position in the sequence
	# seedseq: input sequence withs positions to permute
  # seq: current sequence being built
	# container for adding new sequences to
#	print index # debug
	if seq == []:
		for i in range(len(seedseq)): seq.append( Position() )
	if permute_me( seedseq[index] ):
		for base in bases:
			seq[index] = Position( seedseq[index].chain, seedseq[index].index, seedseq[index].type )
			seq[index].type = threeletter[base]
			if index >= len(seedseq) - 1:
				seqs.append( [ Position( P.chain, P.index, P.type ) for P in seq ] ) # hard copy
#				print [ str(P) for P in seqs[-1] ] # debug
			else: permute_positions( index+1, seedseq, seq, seqs )
	else:
		seq[index] = Position( seedseq[index].chain, seedseq[index].index, seedseq[index].type )
		if index >= len(seedseq) - 1:
			seqs.append( [ Position( P.chain, P.index, P.type ) for P in seq ] ) # hard copy
#			print [ str(P) for P in seqs[-1] ] # debug
		else: permute_positions( index+1, seedseq, seq, seqs )

# think of a way around this code duplication?
def permute_positions_excluding_relatives_of( relseq, index, seedseq, seq, seqs ):
	if seq == []:
		for i in range(len(seedseq)): seq.append( Position() )
	if permute_me( seedseq[index] ):
		for base in bases:
			if threeletter[base] == relseq[index].type: continue
			seq[index] = Position( seedseq[index].chain, seedseq[index].index, seedseq[index].type )
			seq[index].type = threeletter[base]
			if index >= len(seedseq) - 1:
				seqs.append( [ Position( P.chain, P.index, P.type ) for P in seq ] ) # hard copy
#				print [ str(P) for P in seqs[-1] ] # debug
			else: permute_positions_excluding_relatives_of( relseq, index+1, seedseq, seq, seqs )
	else:
		seq[index] = Position( seedseq[index].chain, seedseq[index].index, seedseq[index].type )
		if index >= len(seedseq) - 1:
			seqs.append( [ Position( P.chain, P.index, P.type ) for P in seq ] ) # hard copy
#			print [ str(P) for P in seqs[-1] ] # debug
		else: permute_positions_excluding_relatives_of( relseq, index+1, seedseq, seq, seqs )

def permute_me(P):
	if P.type.lower() in threeletter.keys(): return 0
	if P.type in threeletter.values(): return 0
	return 1

def SNPs(codon):
	codon = codon.upper()
	SNPs = []
	for i in range(len(codon)):
		for base in BASES:
			if codon[i] == base: continue
			# list workaround for 'immutable string' silliness
			SNP = [n for n in codon]
			SNP[i] = base
			SNP = string.join(SNP,'')
			SNPs.append(SNP)
	return SNPs

class Residue:
	def __init__(self,type=''):
		self.type = type
		self.coord = {}
		self.behavior = []
		self.mutable = False
		self.chi_deltas = []

	def __str__(self):
		return self.type

	def add_atom(self,line):
		atom = line[12:16].strip()
#		exclude hydrogens and backbone atoms
		if skip_atom(atom): return
		x = line[30:38].strip(); y = line[38:46].strip(); z = line[46:54].strip()
#		print atom, x, y, z
		self.coord[atom] = ( float(x), float(y), float(z) )
#		self.print_atom(atom)

	def print_atom(self,atom):
		print '%3s %4s %7.3f %7.3f %7.3f' % ( self.type, atom, self.coord[atom][0],
		                                  self.coord[atom][1], self.coord[atom][2] )

class Chain:
	def __init__(self):
		self.id = '-'
		self.residues = {}
	def add(self,chainpos,residue=Residue()):
		if self.id == '-': self.id = chainpos.chain
		elif chainpos.chain != self.id:
			print 'error, attempting to combine chains!'
		self.residues[chainpos] = residue
	def res(self,chainpos):
		if self.residues.has_key(chainpos):
			return self.residues[chainpos]
		else:
			print 'error: no %s in chain!' %chainpos
			return None
	def seq(self):
		seq = []
		for cp in sorted(self.residues):
			seq.append( self.residues[cp].type )
		return string.join( seq, '' )
	def strdefs(self):
		out = []
		for cp in sorted(self.residues):
			out.append( '%s.%i.%s' %(cp.chain, cp.pos, self.residues[cp]) )
		return string.join(out, '\n') + '\n'
	def __str__(self):
		out = 'Chain %s: %s to %s\n' %(self.id, self.start(), self.end())
		return out + self.seq()
	def fasta(self,prefix=''):
		return '>%s%s_%s-%s\n%s\n' %(prefix,self.id, self.start(), self.end(), self.seq())
	def start(self):
		return min(self.residues.keys()).pos
	def end(self):
		return max(self.residues.keys()).pos
	def __cmp__(self,other):
		if self.id == None or other.id == None: return 0
		if self.id < other.id: return -1
		if self.id > other.id: return 1
		return 0

class ChainPos:
	def __init__(self,chain,pos):
		self.chain = chain
		self.pos = int(pos)
	def __str__(self):
		return '%s.%i' %(self.chain,self.pos)
	def __lt__(self,other):
		if self.chain < other.chain: return True
		if self.chain > other.chain: return False
		if self.pos < other.pos: return True
		if self.pos > other.pos: return False
		return False # strict weak ordering
	def __cmp__(self,other):
		if self.chain < other.chain: return -1
		if self.chain > other.chain: return 1
		if self.pos < other.pos: return -1
		if self.pos > other.pos: return 1
		return 0
	def __hash__(self):
		return hash( (self.chain,self.pos) )

def dist3D( xyz1, xyz2 ):
	dx = xyz2[0]-xyz1[0]; dy = xyz2[1]-xyz1[1]; dz = xyz2[2]-xyz1[2];
	return sqrt( dx*dx + dy*dy + dz*dz )

def residue_rmsd(res1,res2):
	l1 = len(res1.coord)
	if len(res2.coord) != l1:
		print 'Warning, rmsd calculation for residues of different length!'
	if l1 == 0: return 0
	sum = 0.0
	for atom in res1.coord:
#		res1.print_atom(atom); res2.print_atom(atom)
		if not res1.coord.has_key( atom ): continue;
		if not res2.coord.has_key( atom ): continue;
		dis = dist3D( res1.coord[atom], res2.coord[atom] )
		sum += dis*dis
	return sqrt(float(sum)/float(l1))

def make_resfile( pdb, resname = '', restag = 'AUTOM' ):
	if not os.path.exists(pdb):
		print '%s not found' %pdb; return
	root = pdb.split('.')[0]
	if resname == '': resname = '%s.res' %root
	resfile = open( resname, 'w' )
	resfile.write(' start\n') # Rosetta resfiles need this
	resj = 0; lastres = -1
	for line in open(pdb,'r'):
		if not line.startswith('ATOM'): continue
		resn = line[17:20]
		chain = line[21:22]
		resi = int( line[22:26].strip() )
		tag = restag
		if resn in nucleic_acids: tag = 'NATVA'
		if resi != lastres:
			resj += 1; lastres = resi
			resfile.write('%2s%5i%5i%6s%4s%2s\n' % \
			              (chain,resj,resi,tag,'#',oneletter[resn]) )
	return resname

def make_resfile_new( pdb, resname = '', global_defaults = ['AUTO'], default = '' ):
	if not os.path.exists(pdb):
		print '%s not found' %pdb; return
	root = pdb.split('.')[0]
	if resname == '': resname = '%s.resfile' %root
	resfile = open( resname, 'w' )
	resfile.write('# generated from %s\n' %pdb )
	resfile.write( '%s\n' % string.join( global_defaults, ' ' ) )
	resfile.write('start\n')
	resj = 0; lastres = -1
	for line in open(pdb,'r'):
		if not line.startswith('ATOM'): continue
		resn = line[17:20]
		chain = line[21:22]
		resi = int( line[22:26].strip() )
		tag = default
		if resn in nucleic_acids: tag = 'NATRO'
		if resi != lastres:
			resj += 1; lastres = resi
			thisoneletter = ''
			try: thisoneletter = oneletter[resn]
			except: pass
			resfile.write('%i %s %s # %s\n' %( resi, chain, tag, thisoneletter ) )

def children( files, prefix ):
	return sorted( [ path for path in files if path[:len(prefix)+1] == prefix+'_' ] )

def get_prefix( path, num_keys ):
	return string.join( path.split('_')[:num_keys], '_' )

def chi_deltas_same( a, b ):
	# compares list of tuples with indices 0 : resid, 1-4 : chi deltas
	# order matters: expects that 'uniform' rosetta PDB files were read
	length = len(a)
	if len(b) < length: length = len(b)
	for i in range(length):
		if a[i][0] != b[i][0]: return False
		j = 1
		while j < len(a[i]):
			d = float(a[i][j]) - float(b[i][j])
			if d < -1 or d > 1: return False # 1-degree cushion
			j += 1
	return True

def zgrep( string, path ):
	return commands.getoutput( "zgrep %s %s" % ( string, path ) )

def compseq(seq):
	return string.join( [ compbase(base) for base in seq ], '' )

def rvs_comp_str(seq):
	return string.join([ compbase(base) for base in reversed(seq) ], '' )

def rvs_comp(seq):
	return [ compbase(base) for base in reversed(seq) ]

def translate_simple(seq):
	pseq = []
	re_codon = re.compile('([acgtnACGTN]{3})')
	for m in re_codon.finditer(seq):
		pseq.append( translate_codon( m.groups()[0] ) )
	return string.join(pseq,'')

def translate(seq):
	re_codon = re.compile( '[a-zA-Z]{3}' )
	frames = ([],[],[])
	i = 0
	while i < len(seq)-2:
		for j in range(3):
			if i+j > len(seq)-3: continue
			codon = seq[i+j:i+j+3]
			if re_codon.match(codon): frames[j].append( translate_codon( codon.lower() ) )
			else: frames[j].append( '-' )
		i += 3
	return frames

# homemade decimal to hexidecimal
def dec2hex( dec ):
	hex = [ str(i) for i in range(10) ]
	hex.extend( ['a','b','c','d','e','f'] )
	places = []
	baserecurse( dec, places, 16 )
	return string.join( [ hex[place] for place in reversed(places) ], '' )

# should work for any base 'n'
def baserecurse( dec, places, base ):
	fact = dec/base
	nextplace = 0
	while nextplace < fact:
		nextplace += 1
		dec -= base
	places.append( dec )
	if nextplace > base-1: hexrecurse( nextplace, places )
	else: places.append( nextplace )

def diff_fwd( str1, str2 ):
	if ( len(str1) - len(str2) ) != 0: return -1
	count = 0
	for i in range( len( str1 ) ):
		if str1[i] != str2[i]: count += 1
	return count

def diff_rvs( str1, str2 ):
	return diff_fwd( str1, rvs_comp_str(str2) )

def mark_mismatches( query, ref ):
	marked = ''
	for i in range( len(query) ):
		# case-insensitive w/ respect to the reference
		if query[i] == string.lower(ref[i]) or query[i] == string.upper(ref[i]):
			marked += string.upper(query[i])
		else: marked += string.lower(query[i])
	return marked

def mean_sd( list ):
	l = len(list)
	if l == 1: return list[0], 0.0
	mean = 0.0
	for i in list: mean += i
	mean = mean/l
	sd = 0.0
	for i in list:
		sd += (i-mean)*(i-mean)
	return mean, sqrt(sd/(l-1))

if __name__ == "__main__":
	print 'Python library--no standalone implementation'
