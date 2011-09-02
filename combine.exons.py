#!/usr/bin/env python
import os,re,string,sys

'''
chr_1	CEG	exon	433077	433745	.	+	.	exonNumber 1; transcriptId 260789
chr_1	CEG	exon	433794	434729	.	+	.	exonNumber 2; transcriptId 260789
chr_1	JGI	exon	300	1153	.	-	.	name "fgenesh1_pg.C_chr_1000001"; transcriptId 867
chr_1	JGI	CDS	300	1153	.	-	0	name "fgenesh1_pg.C_chr_1000001"; proteinId 867; exonNumber 10
'''
re_gff = re.compile('(\w+)\s+\w+\s+exon\s+([0-9]+)\s+([0-9]+)\s+.+\s+([+-])\s+.+transcriptId (\w+)')

'''
bd.1806	bd_10x65	-	13467	13902
bd.1807	bd_10x65	+	15901	16371
bd.1807	bd_10x65	+	16804	16971
bd.1807	bd_10x65	+	17108	17988
'''
re_exon = re.compile('([a-zA-Z0-9-_.]+)\s+(\w+)\s+([+-])\s+([0-9]+)\s+([0-9]+)')

sep = '\t'

class Segment:
	sep = sep
	def __init__(self,seq='',strand='',start=0,end=0):
		self.seq = seq
		self.strand = strand
		self.start = start
		self.end = end
	def header(self):
		return string.join( ['seq','strand','start','end'], self.sep )
	def __str__(self):
		return string.join( [self.seq,self.strand,str(self.start),str(self.end)], self.sep )
	def __cmp__(self,other):
		if self.seq < other.seq: return -1
		if self.seq > other.seq: return 1
		if self.start < other.start: return -1
		if self.start > other.start: return 1
		if self.strand == '+' and other.strand == '-': return -1
		if self.strand == '-' and other.strand == '+': return 1
		if self.end < other.end: return -1
		if self.end > other.end: return 1
		return 0

merged = {}

file = sys.argv[1]
#for m in re_gff.finditer(open(file).read()):
#	seq,start,end,strand,id = m.groups()
for m in re_exon.finditer(open(file).read()):
	id,seq,strand,start,end = m.groups()
	start=int(start)
	end=int(end)
	if not merged.has_key(seq): merged[seq] = {}
	if not merged[seq].has_key(id): merged[seq][id] = Segment(seq,strand,start,end)
	else:
		seg = merged[seq][id]
		if start < seg.start: merged[seq][id].start = start
		if end > seg.end: merged[seq][id].end = end

def byval(a,b):
	if a[1] < b[1]: return -1
	if a[1] > b[1]: return 1
	return 0

print '%s%s%s' %('id',sep,Segment().header())
for seq in sorted(merged):
	for id,seg in sorted( merged[seq].items(), byval ):
		print '%s%s%s' %(id,sep,str(seg))
