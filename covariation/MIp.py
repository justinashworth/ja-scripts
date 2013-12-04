#!/usr/bin/env python
import os,re,string,sys,math

class ColPair:
	sep = '\t'
	def __init__(self,a,b,mi=None,mip=None):
		self.a = a
		self.b = b
		self.mi = mi
		self.mip = mip

	def __str__(self):
		# add 1 to alignment indices, since Python is zero-indexed
		out = '%i%s%i%s%.3g' %(self.a+1,self.sep,self.b+1,self.sep,self.mi)
		if self.mip: out = out + '%s%.4g' %(self.sep,self.mip)
		return out

# read alignment, store as columns
ncol = 0
cols = []
for line in open(sys.argv[1]):
	line = line.strip()
	print line

	if ncol==0:
		ncol = len(line)
		for i in range(len(line)): cols.append([])

	elif not len(line)==ncol:
		print 'Error: line %s: non-conforming length' %line
		sys.exit()

	for i in range(len(line)):
		cols[i].append(line[i])

def px(X):
	ps = {}
	for x in set(X):
		ps[x] = float(len(filter(lambda z: z==x, X)))/len(X)
	return(ps)

def MI(A,B):
	l = len(A)
	pa = px(A); pb = px(B)
	mi = 0
	for a in set(A):
		for b in set(B):
			ab = 0
			for i in range(l):
				if A[i]==a and B[i]==b: ab=ab+1
			pab = float(ab)/l
			#print pab,pa[a],pb[b]
			if pab==0: continue
			mi = mi + pab * math.log(pab/pa[a]/pb[b])
	return mi

################################################################################
#### START ####
################################################################################

pairs = {}

# compute mutual information for all pairs of columns in the alignment (as in Fodor and Aldrich 2004)

for i in range(ncol):
	for j in range(i,ncol):
		if i==j: continue
		#print i,j
		# skip gapped columns (as done by Fodor and Aldrich)
		if '-' in cols[i]: continue
		if '-' in cols[j]: continue
		mi = MI(cols[i],cols[j])
		#print mi
		pairs[(i,j)] = ColPair(i,j,mi=mi)

# next: MIp: background-corrected mutual information (Dunn et al. 2008)

# overall average MI
MIav = 0
for pair in pairs.values(): MIav = MIav + pair.mi
MIav = MIav / len(pairs)
print 'MIav',MIav

# compute average of all MIs for any position/column 'a'
MIax = {}
for i in range(ncol):
	miax = 0
	count = 0
	for j in range(ncol):
		if i==j: continue
		if pairs.has_key( (i,j) ):
			miax = miax + pairs[(i,j)].mi
			count = count + 1
		if pairs.has_key( (j,i) ):
			miax = miax + pairs[(j,i)].mi
			count = count + 1
	if count > 0: MIax[i] = miax / count
	else: MIax[i] = 0

# 'MIp'
for pair in pairs.values():
	# = MI - 'average product correction (APC)'
	pair.mip = pair.mi - MIax[pair.a] * MIax[pair.b] / MIav

outf = open('mip.out','w')
outf.write('i\tj\tmi\tmip\n')
for i in range(ncol):
	for j in range(i,ncol):
		if pairs.has_key((i,j)):
			outf.write('%s\n' %pairs[(i,j)])
		if pairs.has_key((j,i)):
			outf.write('%s\n' %pairs[(j,i)])
	outf.write('%s\n' %str(pair))
outf.close()

# wite matrices (such as for plotting heatmaps in R)
matf_mi = open('mi.mat','w')
for i in range(ncol):
	line = []
	for j in range(ncol):
		if i==j:
			line.append( '%.3f' %(0) )
		elif pairs.has_key((i,j)):
			line.append( '%.3f' %pairs[(i,j)].mi )
		elif pairs.has_key((j,i)):
			line.append( '%.3f' %pairs[(j,i)].mi )
	matf_mi.write( '%s\n' %string.join(line, '\t') )

matf_mip = open('mip.mat','w')
for i in range(ncol):
	line = []
	for j in range(ncol):
		if i==j:
			line.append( '%.4f' %(0) )
		elif pairs.has_key((i,j)):
			line.append( '%.4f' %pairs[(i,j)].mip )
		elif pairs.has_key((j,i)):
			line.append( '%.4f' %pairs[(j,i)].mip )
	matf_mip.write( '%s\n' %string.join(line, '\t') )
