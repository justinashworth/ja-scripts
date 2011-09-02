#!/bin/sh
reference=$1
probes=$2
p=$probes.`basename $reference .fasta`
p=`basename $p .fa`
echo "naive probes: `grep -c ">" $probes.fa`" > $p.log
# fine-grained BLAT to see where designed probes actually match
blat $reference $probes.fa $p.blat -tileSize=10 -fine
# filter out BLAT matches that are imperfect, gapped, or partial
filter-blat-matches.py -R $p.blat > $p.filtered.dat
#
blat-probes.R.sh $p.filtered.dat $probes $p
cat $p.log
