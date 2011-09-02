#!/usr/bin/env python
__author__ = 'Justin Ashworth'

import re,sys,string,os
from AshworthUtil import ChainPos,Residue,BASES,bases,get_oneletter,nucleic_acids

class PDBFile:
	def __init__(self):
		self.residues = {}
		self.readcoords = False

	def __str__( self ):
		out = '%s' % ( self.filename )
		reslist = []
		keys = self.residues.keys()
		keys.sort()
		for chainpos in keys:
			res = self.residues[chainpos]
			reslist.append( '%s %s' %(str(chainpos),res.type) )
		out += 'res: ' + string.join( reslist, ', ' )
		return out

	def residue(self,chain,pos):
		if not self.has_residue(chain,pos):
			print 'warning, no chain %s pos %s in pdb' %(chain,pos)
			return
		return self.residues[ChainPos(chain,pos)]

	def has_residue(self,chain,pos):
		return self.residues.has_key(ChainPos(chain,pos))

	def read_pdb( self, filename, readcoords=False ):
		self.filename=filename
		re_base_atom = re.compile('^ATOM.+ CA |^ATOM.+ C1[*\'] ')
		pdb = None
		if filename.split('/')[-1].split('.')[-1] == 'gz':
			pdb = os.popen( 'gunzip -c %s' % filename )
		else: pdb = file( filename, 'r' )

		print 'reading %s' % filename
		curr_chainpos = ChainPos('-',0)
		for line in pdb:
			if line.startswith('ATOM'):
				chainpos = ChainPos( line[21:22], line[22:26].strip() )
				if re_base_atom.search(line): # usually the second atom, after N
					curr_chainpos = chainpos
					self.residues[chainpos] = Residue( line[17:20] )
				if readcoords and chainpos == curr_chainpos:
						self.residues[chainpos].add_atom(line)

	def fasta( self ):
		os = ['> %s\n' %self.filename]
		# note: self.residues is indexed by instances of ChainPos
		keys = self.residues.keys()
		keys.sort() # ChainPos knows how to sort itself
		for chainpos in keys:
			res = self.residues[chainpos]
			# skip DNA
			if res.type in nucleic_acids: continue
			letter = get_oneletter(res.type)
			os.append(letter)
		os.append('\n')
		return string.join( os, '' )

	def dna(self):
		os = ['> %s\n' %self.filename]
		keys = self.residues.keys()
		keys.sort()
		for chainpos in keys:
			res = self.residues[chainpos]
			if not res.type in nucleic_acids: continue
			letter = get_oneletter(res.type)
			os.append(letter)
		os.append('\n')
		return string.join( os, '' )

### end of PDBFile class

class RosettaPDBFile(PDBFile):
	def __init__( self, input='', readcoords=False, match_list = ['energy','binding','base-spec'] ):
		PDBFile.__init__(self)
		self.filename = input
		self.E = {}

		self.read_pdb( input, readcoords, match_list )

		for term in ['energy','binding','base-spec']:
			if not self.E.has_key(term): self.E[term] = 0.0

	def __str__( self ):
		out = '%s' % ( self.filename )
		for val in self.E.values(): out += '%10.2f' % val
		return out

	def str_with_res( self ):
		out = '%s' % ( self.filename )
		for val in self.E.values():
			out += '%10.2f ' % val
		reslist = []
		keys = self.residues.keys()
		keys.sort()
		for chainpos in keys:
			res = self.residues[chainpos]
			reslist.append( '%s %s' %(str(chainpos),res.type) )
		out += 'res: ' + string.join( reslist, ', ' )
		return out

	def __eq__( self, other ):
		dbkt = self.E['energy'] - other.E['energy']
		dddg = self.E['binding'] - other.E['binding']
		# 0.05 energy unit cushion
		return ( dbkt > -0.05 and dbkt < 0.05 and
		         dddg > -0.05 and dddg < 0.05 )
#		return ( dbkt > -0.05 and dbkt < 0.05 and
#		         dddg > -0.05 and dddg < 0.05 and
#		         chi_deltas_same( self.chi_deltas, other.chi_deltas ) )

	def __cmp__( self, other ):
		d = self.E['energy'] - other.E['energy']
		if d < 0: return -1
		elif d > 0: return 1
		return 0

	# list of mutable residues (those residues that were allowed to mutate)
	def mutable( self ):
		mutable = []
		for chainpos,residue in self.residues.items():
			for type in ['DesignRes','MutatedRes','MutatedDNA']:
				if type in residue.behavior:
					mutable.append(chainpos)
					break
		return mutable

	# each of these is their own function for simplicity
	def line_to_energy(self,line):
		for field in line.split():
			if re.match( '^[0-9.-]+$', field ): self.E['energy'] = float(field)
	def line_to_binding(self,line): self.E['binding'] = float(line.strip().split(':')[-1])
	def line_to_Fitness(self,line): self.E['Fitness'] = float(line.strip().split(':')[-1])
	def line_to_base_spec(self,line): self.E['base-spec'] = float(line.split()[3])
	def line_to_chi_offset(self,line):
		return

	def line_to_molten(self,line):
		comma = line.strip().split(',')
		type = comma[0].split()[-1]
		#print type
		for i in comma[1:]:
			i = i.strip().split()
			if len(i) < 2: continue
			chainpos = ChainPos(i[1],i[0])
			if self.residues.has_key(chainpos): self.residues[chainpos].behavior.append( type )
			else: print 'Warning, %s residue %s not found in pdb!' % ( type, chainpos )

	def read_pdb( self, filename, readcoords, match_list ):

		# make regex for special case
		re_base_atom = re.compile('^ATOM.+ CA |^ATOM.+ C1[*\'] ')

		match_fxns = {
			'molten' : ( re.compile( 'PackingRes|DesignRes|MutatedRes|MovedRes|MutatedDNA|MovedDNA|DNAdesRes' ), self.line_to_molten ),
			'energy' : ( re.compile( 'bk_tot|total_score' ), self.line_to_energy ),
			'binding' : ( re.compile( 'Binding' ), self.line_to_binding ),
			'Fitness' : ( re.compile( 'Fitness' ), self.line_to_Fitness ),
			'base-spec' : ( re.compile( 'base-spec' ), self.line_to_base_spec ),
			'chi_offset' : ( re.compile( 'chi_offset' ), self.line_to_chi_offset ),
		}

		pdb = None
		if filename.split('/')[-1].split('.')[-1] == 'gz':
			pdb = os.popen( 'gunzip -c %s' % filename )
		else: pdb = file( filename, 'r' )

		print 'reading %s' % filename
		curr_chainpos = ChainPos('-',0)
		molten_lines = []
		for line in pdb:
			if line.startswith('ATOM'):
				chainpos = ChainPos( line[21:22], line[22:26].strip() )
				if re_base_atom.search(line):
					curr_chainpos = chainpos
					self.residues[chainpos] = Residue( line[17:20] )
				if readcoords and chainpos == curr_chainpos:
						self.residues[chainpos].add_atom(line)
			else:
				for match in match_list:
					if not match_fxns.has_key( match ): continue
					if match_fxns[match][0].search(line):
						if match == 'molten': molten_lines.append(line); break
						match_fxns[match][1](line)

		for line in molten_lines:
			self.line_to_molten(line)

	def set_specificity( self, value ): self.specificity = value
### end of RosettaPDBFile class

if __name__ == '__main__':
	print 'Python class library--no standalone implementation'
