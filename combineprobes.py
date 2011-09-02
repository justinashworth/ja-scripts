#!/usr/bin/env python
import os,re,string,sys
from optparse import OptionParser
p=OptionParser()
p.add_option('-m','--mapping')
p.add_option('-d','--data')
p.add_option('-s','--sep',default='\t')
opt,args=p.parse_args()

data={}
colnames=[]
for l in open(opt.data):
	words=l.strip().split()
	if l.startswith('probe'): 
		colnames=words
		continue
	probe=words[0]
	data[probe]={}
	for i in range(len(words)-1):
		data[probe][colnames[i+1]]=words[i+1]

mapping={}
for l in open(opt.mapping):
	if l.startswith('probe'): continue
	words=l.strip().split()
	probe=words[0]
	if not probe in data.keys(): continue
	id=words[1]
	if not mapping.has_key(id): mapping[id]=[]
	mapping[id].append(probe)

def asint(a,b):
	return int(a) - int(b)

outfname='%s.geneaveraged'%opt.data
if os.path.exists(outfname): os.remove(outfname)
outf=open(outfname,'w')
outf.write('%s%s%s\n' %('id',opt.sep,string.join(colnames[1:],opt.sep)))
ids=mapping.keys()
ids.sort(asint)
for id in ids:
	#print id
	probes=mapping[id]
	tempcols={}
	for col in colnames[1:]:
#		print col
		tempcols[col]=[0,0]
		for probe in probes:
#			print probe
			if not data.has_key(probe): continue
			tempcols[col][0]+=float(data[probe][col])
			tempcols[col][1]+=1
	outf.write('%s%s%s\n' %(id,opt.sep,string.join(['%.3f' %float(tempcols[c][0]/tempcols[c][1]) for c in colnames[1:]],opt.sep)))	
