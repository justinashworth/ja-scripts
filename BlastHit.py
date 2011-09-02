#!/usr/bin/env python
import os,re,string,sys
from Fasta import *
from optparse import OptionParser

__author__ = "Justin Ashworth"

'''
DESCRIPTION
	Basic parser for standard BLAST output. Work in progress. Originally written mainly to convert annoying multiline (wrapped) alignment format to something remotely usable for downstream applications.
	Note: no mechanism yet to maintain master alignment (needs to take gap differences between different hits into account)
'''

'''
Example of format to be parsed:
'''
'''
> scaffold3858 36.8
Length=1460

 Score = 37.4 bits (85),  Expect = 3e-04, Method: Compositional matrix adjust.
  Identities = 21/54 (38%), Positives = 36/54 (66%), Gaps = 3/54 (5%)
 Frame = -2

 Query  34   VDEAENGKEALAKLKQQKYDLVLLDWNMPEMDGITLLQEIKKDPDLKNIPVIMV  87
             V +++NG +A  K+++   DL+LLD  +P ++GI L ++IK   D   +P+IMV
 Sbjct  610  VTKSDNGIDAWNKIEKNLPDLILLDIMIPGINGIDLCKKIKSQYD---VPIIMV  458

> C91383 10.0
 Length=109

  Score = 34.7 bits (78),  Expect = 0.002, Method: Composition-based stats.
 Identities = 13/36 (36%), Positives = 24/36 (66%), Gaps = 0/36 (0%)
  Frame = +1

Query  49   QQKYDLVLLDWNMPEMDGITLLQEIKKDPDLKNIPV  84
            +++  L++ D  MP MDG  +  EI+K+P  ++IP+
Sbjct  1    RERPALIVADILMPRMDGFAMAHEIRKNPQTRSIPI  108
'''

#class BLASTSeq(Fasta):
#	def __init__(self,seq='',name='',type='',start=0,end=0):
#		Fasta.__init__(self,seq,name,type)
#		self.start = int(start)
#		self.end = int(end)
#	def append(self,seq,start

class BLASTHit:
	re_hit = re.compile('>(.*?)\n.*?(?=>)',re.DOTALL)
#	re_hit = re.compile('>(.*?)\nLength=([0-9]+).*?Expect = ([0-9.-e]+).*?Query',re.DOTALL)
	re_query = re.compile('Query +([0-9]+) +([A-Za-z-]+) +([0-9]+)')
	re_subject = re.compile('Sbjct +([0-9]+) +([A-Za-z-]+) +([0-9]+)')
#	re_subject = re.compile()

	def __init__(self,match):
		matched = match.group()
		self.q_start=0
		self.q_end=0
		self.q_seq=''
		self.s_start=0
		self.s_end=0
		self.s_seq=''

		subject = match.groups()[0]
		self.subject = re.sub('\n',' ',subject.strip())

		for s in self.re_query.finditer(matched):
			start,seq,end = s.groups()
			start=int(start); end=int(end)
			if self.q_start == 0: self.q_start = start
#			if end > self.q_end: self.q_end = end
			self.q_seq += seq

		for s in self.re_subject.finditer(matched):
			start,seq,end = s.groups()
			start=int(start); end=int(end)
			if self.s_start == 0: self.s_start = start
#			if end > self.s_end: self.s_end = end
			self.s_seq += seq

		if len(self.q_seq) != len(self.s_seq):
			sys.stderr.write('warning: sequence length mismatch (%s)\n' %self.subject)

#	if self.q_end-self.q_start != len(self.q_seq):
#		sys.stderr.write('warning: query sequence length didn\'t match annotated numbering\n')
#	if self.s_end-self.s_start != len(self.s_seq):
#		sys.stderr.write('warning: subject sequence length didn\'t match annotated numbering\n')

	def expand_query(self,query_length=0):
		leading = string.join( [ '-' for i in range(self.q_start-1) ], '' )
		self.q_seq = leading + self.q_seq
		self.s_seq = leading + self.s_seq
		self.q_start = 1
		self.s_start = 1
#		lagging = []
#		cursor = self.q_end
#		while cursor <= query_length:
#			lagging.append('-')
#			cursor += 1
#		lagging = string.join( lagging, '' )
#		self.q_seq = self.q_seq + lagging
#		self.s_seq = self.s_seq + lagging
#		self.q_end += len(lagging)
#		self.s_end += len(lagging)

	def __str__(self):
		out = []
		out.append( '%s\t%i\t%s' %(self.subject,self.s_start,self.s_seq) )
		out.append( 'Query\t%i\t%s' %(self.q_start,self.q_seq) )
		return string.join(out,'\n')

	def subject_fa(self):
		out = []
		out.append( '>%s_%i\n%s' %(self.subject,self.s_start,self.s_seq) )
		return string.join(out,'\n')

class BLASTParse:
	re_query = re.compile('Query=(.*?)Length=(.*?)\n',re.DOTALL)
	def __init__(self):
		self.hits = []

	def run(self,opt,args):
		query_length = 0
		for hitfile in args:
			if not os.path.exists(hitfile):
				sys.stderr.write('warning, hitfile %s not found\n' %hitfile)
				continue
			qmatch = self.re_query.search(open(hitfile).read())
			if not qmatch:
				sys.stderr.write('warning, no query id found in %s\n' %hitfile)
				continue
			query,query_length = qmatch.groups()
			query = re.sub('\n',' ',query.strip())
			query_length = int(query_length)
			sys.stderr.write('Query: %s (length: %s)\n' %(query,query_length))
			for match in BLASTHit.re_hit.finditer(open(hitfile).read()):
				self.hits.append( BLASTHit(match) )
		self.expand_query(query_length)

	def expand_query(self,query_length):
		for h in self.hits:
			h.expand_query(query_length)

	def __str__(self):
		out = []
#		for h in self.hits: out.append( str(h) )
		for h in self.hits: out.append( h.subject_fa() )
		return string.join(out,'\n')

if __name__ == "__main__":
	p=OptionParser()
	opt,args=p.parse_args()
	app = BLASTParse()
	app.run(opt,args)
	print str(app)
