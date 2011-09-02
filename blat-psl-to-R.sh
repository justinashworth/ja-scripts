#!/bin/bash
outf=`basename $1 .psl`.tsv
# header for genes
#echo 'match	mismatch	repmatch	Ns	Qgapcount	Qgapbases	Tgapcount	Tgapbases	strand	Qid	Qname	Qsize	Qstart	Qend	Tname	Tsize	Tstart	Tend' > $outf
# header for probes (no 'Qid' field)
echo 'match	mismatch	repmatch	Ns	Qgapcount	Qgapbases	Tgapcount	Tgapbases	strand	Qname	Qsize	Qstart	Qend	Tname	Tsize	Tstart	Tend' > $outf
#cat $1 | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' | tail -n +6 | sed "s/jgi|Thaps3|//" | tr '|' '\t' | sed "s/	[Cc]hr_\([0-9]*[a-z]\)_[0-9]*/	Chr\1/" | sed "s/	chr_/	Chr/" >> $outf
cat $1 | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17}' | tail -n +6 >> $outf
