#!/usr/bin/env python
import os,re,string,sys
from optparse import OptionParser
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import Element

'''
DESCRIPTION
	XML NCBI gene summary file to table column format
	MEMORY-INEFFICIENT: parser as implemented loads entire file into huge memory footprint. There ought to be a memory-efficient iteration method to this, but I haven't seen it yet.
	Beware HUGE XML files even for reasonable nubmers (~10000) of genes. XML files over 100MB may need to be split up and processed separately.

	(NOTE: I have since started using asn2flat utility from NCBI c++ toolkit to accomplish basically the same thing, though this script still provides convenient XML parsing for custom tabulations of various info)
'''

class Gene:
	sep = '\t'
	def __init__(self,name='',tag='',gid='',genome='',start=0,end=0,strand='',source='',protein=''):
		self.name=name
		self.tag=tag
		self.gid=gid
		self.genome=genome
		self.start=int(start)
		self.end=int(end)
		self.strand=strand
		self.source=source
		self.protein=protein

	def parse(self,gene):
		# 'gene' is an Element
			geneid = gene.findtext('Entrezgene_track-info/Gene-track/Gene-track_geneid')
			if geneid: self.gid = geneid
			genome = gene.find('Entrezgene_source/BioSource/BioSource_genome')
			if genome != None:
				value = genome.get('value')
				if value != 'chromosome': self.genome = value
				else: self.genome = 'chr_' + gene.findtext('Entrezgene_source/BioSource/BioSource_subtype/SubSource/SubSource_name')
			else:
				# this seems to apply to not-fully-assembled Thaps reference contigs, e.g. 'chr_16a'
				accession = gene.findtext('Entrezgene_locus/Gene-commentary/Gene-commentary_accession')
				if accession: self.genome = accession
			source = gene.findtext('Entrezgene_source/BioSource/BioSource_org/Org-ref/Org-ref_taxname')
			if source: self.source = source
			name = gene.findtext('Entrezgene_gene/Gene-ref/Gene-ref_locus')
			if name: self.name = name
			tag = gene.findtext('Entrezgene_gene/Gene-ref/Gene-ref_locus-tag')
			if tag: self.tag = tag
			protein = gene.findtext('Entrezgene_prot/Prot-ref/Prot-ref_name/Prot-ref_name_E')
			if protein: self.protein = protein
			locus = gene.find('Entrezgene_locus/Gene-commentary/Gene-commentary_seqs/Seq-loc/Seq-loc_int/Seq-interval')
			if locus:
				self.start = int( locus.findtext('Seq-interval_from') )
				self.end = int( locus.findtext('Seq-interval_to') )
				self.strand = locus.find('Seq-interval_strand/Na-strand').get('value')

	def header(self):
		return string.join(['name','tag','gid','genome','start','end','strand','source','protein'],self.sep)
	def __str__(self):
		return string.join([self.name,self.tag,self.gid,self.genome,str(self.start),str(self.end),self.strand,self.source,self.protein],self.sep)

class App:

	def run(self,opt,files):
		print Gene().header()
		for f in files:
			tree = ET.parse(f)
			iterL = tree.findall('Entrezgene')
			for gene in iterL:
				g = Gene()
				g.parse(gene)
				print str(g)

if __name__ == "__main__":
	p=OptionParser()
	opt,args=p.parse_args()
	app = App()
	app.run(opt,args)
