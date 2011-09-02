#!/usr/bin/env python
import sys,string
from sys import stdout
from multiprocessing import Pool, cpu_count, Manager
from optparse import OptionParser

op=OptionParser()
op.add_option('-l','--lengths',default='6,8')
opt,args=op.parse_args()

# Generates all possible sequences with letters to the length of depth.
def allKmers(depth,letters=['A','C','G','T'],seqs=[''],curdepth=0):
    newseqs = []
    for seq in seqs:
        for letter in letters:
            newseqs.append(seq + letter)
    if depth > curdepth:
        return(allKmers(depth,letters,newseqs,curdepth + 1))
    else:
        return(seqs)

# Gets overlapping number of string instances
def overcount(seq, pattern):
    ocu=0
    x=0
    while 1:
        try:
            i = seq.index(pattern,x)
        except ValueError:
            break
        ocu+=1
        x=i+1
    stdout.write('.')
    stdout.flush()
    return ocu

def ocPool(pattern):
    return overcount(p3utrSeqs[0],pattern)

class KMerHashCount(dict):

	# class to walk through sequences, tallying kmers
	sep = '  '
	def __init__(self,kmers=[],seqs=[],lengths=[6],count_unexpected=True):
		dict.__init__(self)
		self.init(kmers)
		self.count_kmers(seqs,lengths,count_unexpected)

	def init(self,kmers=[]):
		# optional: this basically just makes sure there will be a 'zero' count for expected kmers that didn't occur
		for kmer in kmers:
			self[kmer] = 0

	def count_kmers(self,seqs=[],lengths=[6],count_unexpected=True):
		sys.stderr.write('. = 1 sequence\n')
		for seq in seqs:
#			sys.stderr.write('seq: %s...\n' %seq[:10])
			sys.stderr.write('.')
			l = len(seq)
			i = 0
			while i < l:
				for length in lengths:
					if i+length >= l: continue
					subseq = seq[i:i+length].upper()
					if not self.has_key(subseq):
						sys.stderr.write('\nunexpected k-mer %s encountered...' %subseq)
						if count_unexpected:
							sys.stderr.write('adding to hash')
							self[subseq] = 0
						else:
							sys.stderr.write('skipping')
							continue
					self[subseq] += 1
					# need to hash reverse complement k-mers, or not?
				i += 1
		sys.stderr.write('\n')

	def __str__(self):
		return self.output_std()

	def output_std(self):
		out = []
		keys = self.keys()
		keys.sort()
		for kmer in keys:
			out.append('%s%s%i' %(kmer,self.sep,self[kmer]))
		return string.join(out,'\n')

###########
## p3utr ##
###########
# 1. Read in all sequnces and concatenate them
mgr = Manager()
print 'Reading in sequences...'
#inFile = open('p3utrSeqs_Homo_sapiens.csv','r')
#filename = 'p3utrSeqs_Homo_sapiens.csv'
filename = args[0]
inFile = open(filename,'r')
#seqs3pUTR = [line.strip().split(',')[1] for line in inFile.readlines() if line.strip().split(',')[1]]
# default iterator for a file object is same as .readlines()
seqs3pUTR = [line.strip().split(',')[1] for line in inFile if line.strip().split(',')[1]]
inFile.close()

#ja (I skipped this part)
import numpy
seqLens = [float(len(i)) for i in seqs3pUTR]
p3utrMedLen = int(numpy.median(seqLens))
for i in range(len(seqs3pUTR)):
    if len(seqs3pUTR[i])>p3utrMedLen:
        seqs3pUTR[i] = (seqs3pUTR[i])[1:p3utrMedLen]
p3utrSeqs = mgr.list(['X'.join(seqs3pUTR)])

lengths = opt.lengths.split(',')
for length in lengths:
	length = int(length)
	print 'Making %imers...\n' %length
	# sanity check
	if length < 1 or length > 20:
		sys.stderr.write('are you sure you really wanted kmers of length %i??? skipping\n' %length)
		continue
	kmers = allKmers(length)
	outfilename = '%s.%i.freq' %(filename,length)
	outFile = open(outfilename,'w')
# could also be multiprocessed as multiple KMerHashCount objects, followed by adding up hashes (or 'global' shared hash across multiple walkers) (TO-DO but didn't seem to be necessary for ~10,000 800bp UTRs)
	counts = KMerHashCount(kmers,seqs3pUTR,[length],False)
	outFile.write(str(counts))
	outFile.close()

## 2. Get all 6mers
#print 'Making 6mers...'
#permKMers_6mer = allKmers(6)
#
## 3. Count the instances of all 6mers
#print 'Counting 6mers...'
##outfilename = 'HS3PV2.6.freq'
#outfilename = '%s.6.freq' %filename
#outFile = open(outfilename,'w')
#
#cpus = cpu_count()
#pool = Pool(processes=cpus)
#counts = pool.map(ocPool,permKMers_6mer)
#outFile.write('\n'.join(['  '.join([permKMers_6mer[i],str(counts[i])]) for i in range(len(counts))]))
#outFile.close()
print 'Done.'

## 4. Get all 8mers
#print 'Making 8mers...'
#permKMers_8mer = allKmers(8)
#
## 3. Count the instances of all 6mers
#print 'Counting 8mers...'
##outfilename = 'HS3PV2.8.freq'
#outfilename = '%s.8.freq' %filename
#outFile = open(outfilename,'w')
#
#cpus = cpu_count()
#pool = Pool(processes=cpus)
#
##counts = pool.map(ocPool,permKMers_8mer)
##outFile.write('\n'.join(['  '.join([permKMers_8mer[i],str(counts[i])]) for i in range(len(counts))]))
#
#outFile.close()
#print 'Done.'
"""
###########
## p5utr ##
###########
# 1. Read in all sequnces and concatenate them
print 'Reading in sequences...'
inFile = open('p5utrSeqs_Homo_sapiens.csv','r')
seqs5pUTR = [line.strip().split(',')[1] for line in inFile.readlines()]
inFile.close()
p5utrSeqs = 'X'.join(seqs5pUTR)

# 2. Get all 6mers
print 'Making 6mers...'
permKMers_6mer = allKmers(6)

# 3. Count the instances of all 6mers
print 'Counting 6mers...'
outFile = open('HS5P.6.freq','w')
tmp = permKMers_6mer.pop(0)
outFile.write(str(tmp)+'  '+str(overcount(p5utrSeqs,tmp)))
cur = 1
for kmer in permKMers_6mer:
    outFile.write('\n'+str(kmer)+'  '+str(overcount(p5utrSeqs,kmer)))
    cur += 1
    if cur%100==0:
        stdout.write(str(cur))
    elif cur%10==0:
        stdout.write('.')
    stdout.flush()
outFile.close()

# 4. Get all 8mers
print 'Making 8mers...'
permKMers_8mer = allKmers(8)

# 3. Count the instances of all 6mers
print 'Counting 8mers...'
outFile = open('HS5P.8.freq','w')
tmp = permKMers_8mer.pop(0)
outFile.write(str(tmp)+'  '+str(overcount(p5utrSeqs,tmp)))
cur = 1
for kmer in permKMers_8mer:
    outFile.write('\n'+str(kmer)+'  '+str(overcount(p5utrSeqs,kmer)))
    cur += 1
    if cur%100==0:
        stdout.write(str(cur))
    elif cur%10==0:
        stdout.write('.')
    stdout.flush()
outFile.close()
print 'Done.'

#########
## cds ##
#########
# 1. Read in all sequnces and concatenate them
print 'Reading in sequences...'
inFile = open('cdsSeqs_Homo_sapiens.csv','r')
seqsCDS = [line.strip().split(',')[1] for line in inFile.readlines()]
inFile.close()
cdsSeqs = 'X'.join(seqsCDS)

# 2. Get all 6mers
print 'Making 6mers...'
permKMers_6mer = allKmers(6)

# 3. Count the instances of all 6mers
print 'Counting 6mers...'
outFile = open('HSCDS.6.freq','w')
tmp = permKMers_6mer.pop(0)
outFile.write(str(tmp)+'  '+str(overcount(cdsSeqs,tmp)))
cur = 1
for kmer in permKMers_6mer:
    outFile.write('\n'+str(kmer)+'  '+str(overcount(cdsSeqs,kmer)))
    cur += 1
    if cur%100==0:
        stdout.write(str(cur))
    elif cur%10==0:
        stdout.write('.')
    stdout.flush()
outFile.close()

# 4. Get all 8mers
print 'Making 8mers...'
permKMers_8mer = allKmers(8)

# 3. Count the instances of all 6mers
print 'Counting 8mers...'
outFile = open('HSCDS.8.freq','w')
tmp = permKMers_8mer.pop(0)
outFile.write(str(tmp)+'  '+str(overcount(cdsSeqs,tmp)))
cur = 1
for kmer in permKMers_8mer:
    outFile.write('\n'+str(kmer)+'  '+str(overcount(cdsSeqs,kmer)))
    cur += 1
    if cur%100==0:
        stdout.write(str(cur))
    elif cur%10==0:
        stdout.write('.')
    stdout.flush()
outFile.close()
print 'Done.'

#################
## cds + p3utr ##
#################
# 1. Read in all sequnces and concatenate them
print 'Reading in sequences...'
inFile = open('cdsSeqs_Homo_sapiens.csv','r')
seqsCDS = [line.strip().split(',')[1] for line in inFile.readlines()]
inFile.close()
inFile = open('p3utrSeqs_Homo_sapiens.csv','r')
seqs3pUTR = [line.strip().split(',')[1] for line in inFile.readlines()]
inFile.close()
seqs = []
for i in range(len(seqsCDS)):
    seqs.append(str(seqsCDS[i])+str(seqs3pUTR[i]))
cdsSeqs = 'X'.join(seqs)

# 2. Get all 6mers
print 'Making 6mers...'
permKMers_6mer = allKmers(6)

# 3. Count the instances of all 6mers
print 'Counting 6mers...'
outFile = open('HSMRG.6.freq','w')
tmp = permKMers_6mer.pop(0)
outFile.write(str(tmp)+'  '+str(overcount(cdsSeqs,tmp)))
cur = 1
for kmer in permKMers_6mer:
    outFile.write('\n'+str(kmer)+'  '+str(overcount(cdsSeqs,kmer)))
    cur += 1
    if cur%100==0:
        stdout.write(str(cur))
    elif cur%10==0:
        stdout.write('.')
    stdout.flush()
outFile.close()

# 4. Get all 8mers
print 'Making 8mers...'
permKMers_8mer = allKmers(8)

# 3. Count the instances of all 6mers
print 'Counting 8mers...'
outFile = open('HSMRG.8.freq','w')
tmp = permKMers_8mer.pop(0)
outFile.write(str(tmp)+'  '+str(overcount(cdsSeqs,tmp)))
cur = 1
for kmer in permKMers_8mer:
    outFile.write('\n'+str(kmer)+'  '+str(overcount(cdsSeqs,kmer)))
    cur += 1
    if cur%100==0:
        stdout.write(str(cur))
    elif cur%10==0:
        stdout.write('.')
    stdout.flush()
outFile.close()
print 'Done.'

#########################
## p5utr + cds + p3utr ##
#########################
# 1. Read in all sequnces and concatenate them
print 'Reading in sequences...'
inFile = open('p5utrSeqs_Homo_sapiens.csv','r')
seqs5pUTR = [line.strip().split(',')[1] for line in inFile.readlines()]
inFile.close()
inFile = open('cdsSeqs_Homo_sapiens.csv','r')
seqsCDS = [line.strip().split(',')[1] for line in inFile.readlines()]
inFile.close()
inFile = open('p3utrSeqs_Homo_sapiens.csv','r')
seqs3pUTR = [line.strip().split(',')[1] for line in inFile.readlines()]
inFile.close()
seqs = []
for i in range(len(seqsCDS)):
    seqs.append(str(seqs5pUTR[i])+str(seqsCDS[i])+str(seqs3pUTR[i]))
cdsSeqs = 'X'.join(seqs)

# 2. Get all 6mers
print 'Making 6mers...'
permKMers_6mer = allKmers(6)

# 3. Count the instances of all 6mers
print 'Counting 6mers...'
outFile = open('HSALL.6.freq','w')
tmp = permKMers_6mer.pop(0)
outFile.write(str(tmp)+'  '+str(overcount(cdsSeqs,tmp)))
cur = 1
for kmer in permKMers_6mer:
    outFile.write('\n'+str(kmer)+'  '+str(overcount(cdsSeqs,kmer)))
    cur += 1
    if cur%100==0:
        stdout.write(str(cur))
    elif cur%10==0:
        stdout.write('.')
    stdout.flush()
outFile.close()

# 4. Get all 8mers
print 'Making 8mers...'
permKMers_8mer = allKmers(8)

# 3. Count the instances of all 6mers
print 'Counting 8mers...'
outFile = open('HSALL.8.freq','w')
tmp = permKMers_8mer.pop(0)
outFile.write(str(tmp)+'  '+str(overcount(cdsSeqs,tmp)))
cur = 1
for kmer in permKMers_8mer:
    outFile.write('\n'+str(kmer)+'  '+str(overcount(cdsSeqs,kmer)))
    cur += 1
    if cur%100==0:
        stdout.write(str(cur))
    elif cur%10==0:
        stdout.write('.')
    stdout.flush()
outFile.close()
print 'Done.'
"""
