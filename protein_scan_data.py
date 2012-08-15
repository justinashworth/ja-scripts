#!/usr/bin/env python

import os,sys,re,string
import logging
from AshworthUtil import get_oneletter
from optparse import OptionParser
import cPickle as pickle
from PositionSubstitution import *

debugoutput = False
debugoutput = True

p = OptionParser()
p.add_option('-p','--pathfile')
p.add_option('-k','--pickle')
opt,args = p.parse_args()

re_substitution_new = re.compile( '.*Scores for mutation to ([a-zA-Z]+) at ([\.A-Z0-9]+)\.([a-zA-Z]+): bound = ([\-.0-9]+) .+ binding = ([\-.0-9]+) .+ specificity.bound = ([\-.0-9]+) .+ specificity.binding = ([\-.0-9]+)' )
re_subs_2012 = re.compile( '.*Scores for sequence state [A-Z]\.[0-9]+.([a-zA-Z]+) at ([A-Z])\.([0-9]+).([a-zA-Z]+): bound = ([\-.0-9]+) .+ binding = ([\-.0-9]+) .+ specificity.bound = ([\-.0-9]+) .+ specificity.binding = ([\-.0-9]+)' )
# Scanning protein positions that interface with DNA position(s) CHAIN.###.TYPE
re_scandna = re.compile('Scanning protein positions that interface with DNA position\(s\) ([A-Z]).([0-9]+).([a-zA-Z0-9]+)')

def sort_chainpos(a,b):
	chaina,posa = a.split('.')
	chainb,posb = b.split('.')
	if chaina < chainb: return -1
	if chaina > chainb: return 1
	if int(posa) < int(posb): return -1
	if int(posa) > int(posb): return 1
	return 0

def parse_outfiles(pathfile):
	complexes = {}

	for line in open(pathfile):
		words = line.strip().split()
		if not len(words) == 2:
			logging.warning('skipping non-conforming line: %s' %line)
			continue
		complex,path = words

		print path
		if not complexes.has_key( complex ): complexes[ complex ] = {}

		scandna = re_scandna.search( open(path).read() )
		if scandna:
			chain,pos,type = re_scandna.search( open(path).read() ).groups()
		else:
			chain,pos,type = ['NA','0','NA']
		DNApos = '%s.%s.%s' %(chain,pos,type)
		print DNApos

		if debugoutput: print complex, DNApos
		if not complexes[ complex ].has_key( DNApos ): complexes[ complex ][ DNApos ] = {}

		for line in open(path):
			line = line.strip()
			#match = re_substitution_new.match(line)
			#if match:
			#	aa, chainpos, nataa, bound, binding, specbound, specbinding = match.groups()
			match = re_subs_2012.match(line)
			if match:
				aa, chain, pos, nataa, bound, binding, specbound, specbinding = match.groups()
				aa = get_oneletter(aa)
				nataa = get_oneletter(nataa)
				if debugoutput:
					print aa, chain, pos, nataa
				chainpos = '%s.%s' %(chain,pos)
				if not complexes[ complex ][ DNApos ].has_key(chainpos):
					complexes[ complex ][ DNApos ][ chainpos ] = Position(nataa)

				complexes[ complex ][ DNApos ][ chainpos ].append( \
					Substitution( aa, bound, binding, specbound, specbinding ) )
	return complexes

def flatfiles(prefix):

	subfile = open('%s.substitutions'%prefix,'w')
	natfile = open('%s.natives'%prefix,'w')

	sep = Substitution.sep
	header = ['complex','dna','prot_chain','prot_position','mut','wt','is_native']
	natheader = ['complex','dna','prot_chain','prot_position','wt','ndestab','ndestab_binding']

	class ScoreType:
		def __init__(self,name,temp=0.):
			self.name = name
			self.temp = temp

	score_types = [
		ScoreType('bound',1),
	#	ScoreType('binding',float(1.0)/0.6), # standard thermodynamics (kcal)
		ScoreType('binding',1.0), # effectively no temperature factor
	#	ScoreType('spec_bound',-10),
		ScoreType('spec_binding',-10)
	]

	for st in score_types:
		header.append( '%s' %st.name )
		header.append( '%s_delta' %st.name )
		if re.search('spec',st.name):
			natheader.append( 'dev_%s' %st.name )
			natheader.append( 'sd_%s' %st.name )
			natheader.append( 'weighted_z_score_%s' %st.name )
		else:
			natheader.append( 'Boltzmann_%s' %st.name )
			natheader.append( 'sd_%s' %st.name )
	natheader.append( 'binding_minus_gly' )
	natheader.append( 'spec_minus_gly' )

	for st in score_types:
		header.append( '%s_vs_best' %st.name )

	subfile.write( string.join(header,sep) +'\n')
	natfile.write( string.join(natheader,sep) +'\n')

	counter = 1
	for complex in sorted( complexes ):
		if debugoutput: print complex
		for dna_pos in sorted( complexes[complex] ):
			if debugoutput: print dna_pos
			for chainpos in sorted( complexes[complex][dna_pos], sort_chainpos ):
				if debugoutput: print chainpos
				chain,pos = chainpos.split('.')

				sublist = complexes[complex][dna_pos][chainpos]
				if not sublist.has_native(): continue
				native = sublist.native()

				# print substitution data line to substitutions file
				# also count number of destabilizing mutations above some energy threshold
				destabilizing = {'bound':2,'binding':1}
				n_destabilizing = {'bound':0,'binding':0}
				for sub in sorted( sublist ):
					dataline = [ complex, dna_pos, chain, pos, sub['aa'], native['aa'], '%i' %(sub['aa'] == native['aa']) ]
					for st in score_types:
						dataline.append( '%.3f' %( sub[st.name] ) )
						ddG = sub[st.name] - native[st.name]
						dataline.append( '%.3f' %( ddG ) )
						if st.name in destabilizing:
							if ddG > destabilizing[st.name]:
								n_destabilizing[st.name] += 1

					for st in score_types:
						dataline.append( '%.3f' %( sub[st.name] - sublist.best(st) ) )
					subfile.write( string.join(dataline,sep) + '\n' )

				# print native optimization data line to native file
				natline = [ complex, dna_pos, chain, pos, native['aa'] ]
				natline.append( '%i' %n_destabilizing['bound'] )
				natline.append( '%i' %n_destabilizing['binding'] )

				for st in score_types:
					mean, sd = sublist.mean_sd(st.name)
					if re.search('spec',st.name):
						dev = 0.; z = 0.
						if sd != 0:
							# deviation
							dev = native[st.name] - mean
							# z-score
							z = dev / sd
							# weight by how specific the position was in the native context
	#						z = z * native[st.name]
							# weight by maximum possible specificity at the position
							z = z * sublist.best(st)
						natline.append( '%.3f' %dev )
						natline.append( '%.3f' %sd )
						natline.append( '%.3f' %z )
					else:
						# Boltzmann: for comparing energy terms
						boltz = sublist.Boltzmann(st)
						natline.append( '%.3f' %boltz )
						natline.append( '%.3f' %sd )
				# report the "classic binding energy per residue" as native - gly
				#natline.append( '%.3f' %( native['binding'] - sublist.aa('GLY')['binding'] ) )
				#natline.append( '%.3f' %( native['spec_binding'] - sublist.aa('GLY')['spec_binding'] ) )
				natline.append( '%.3f' %( native['binding'] - sublist.aa('G')['binding'] ) )
				natline.append( '%.3f' %( native['spec_binding'] - sublist.aa('G')['spec_binding'] ) )
				natfile.write( string.join(natline,sep) + '\n' )

				counter+=1

################################################################################
complexes = None
if opt.pathfile:
	complexes = parse_outfiles(opt.pathfile)
	pickle.dump( complexes, open('position.data.p','wb') )
	flatfiles(opt.pathfile)
elif opt.pickle:
	complexes = pickle.load(open(opt.pickle,'rb'))
	flatfiles(opt.pickle)
