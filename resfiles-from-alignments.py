#!/usr/bin/env python
import os,re,string,sys
from optparse import OptionParser
from PDBFile import *
from Alignments import *
from AshworthUtil import isDNA

'''
DESCRIPTION
Write Rosetta resfiles for 'packer'-based threading by simple amino acid substitution.
'''

class ResfilesFromAlignment:
#	def __init__(self):
	def run(self,opt,args):
		if not opt.pdbtemplate:
			print 'error, must specific pdb template (-p option)'; sys.exit()
		pdbtemplate = opt.pdbtemplate
		if not opt.alignment:
			print 'error, must specific alignment (-a option)'; sys.exit()
		alignmentfile = opt.alignment

		pdb = PDBFile()
		pdb.read_pdb(pdbtemplate)
		print str(pdb)
		alignment = Alignment()
		alignment.run(alignmentfile)
		print str(alignment)

		# alignment suitability check
		if len(alignment.templates) < 1:
			print 'error: no template(s)'; return
		if len(alignment.models) < 1:
			print 'error: no model(s)'; return

		print 'writing resfiles for template %s' %pdbtemplate
		defaultresfile = [
			'AUTO',
			'NATAA',
			'start',
		]
		for t in alignment.templates.values():
			for m in alignment.models.values():
				resfilelines = defaultresfile[:]
				for i in range(len(t.chainseqs)):
					start = t.chainseqs[i].start
					tseq = t.chainseqs[i].sequence
					mseq = m.chainseqs[i]
					if len(tseq) != len(mseq):
						print 'error, model unequal length with template (incompatible):\n%s\n%s\n' %(str(t),str(m)); sys.exit()
					# iterate through template alignment positions, add mutation phrases
					for i in range(len(tseq)):
						pos = start.pos+i
						if not pdb.has_residue(start.chain,pos):
							print 'error, template pdb has no residue position %s.%i' %(start.chain,pos); sys.exit()
						modelaa = mseq[i]
						resfilelines.append('%i %s PIKAA %s' %(pos,start.chain,modelaa))

				keys = pdb.residues.keys()
				keys.sort()
				for chainpos in keys:
					if not isDNA(pdb.residues[chainpos].type): continue
					resfilelines.append('%s %s NATRO' %(chainpos.pos,chainpos.chain))
				fname = '%s.%s.%s' %(pdbtemplate.split('.')[0],t.name,m.name)
				resfile = open('%s.resfile' %fname,'w')
				resfile.write( string.join(resfilelines,'\n') )

if __name__ == "__main__":
	p=OptionParser()
	p.add_option('-a','--alignment',help='alignment file')
	p.add_option('-p','--pdbtemplate',help='template pdb file (must be compatible with alignment template(s))')
	opt,args=p.parse_args()
	app = ResfilesFromAlignment()
	app.run(opt,args)
