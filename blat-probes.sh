#!/bin/sh
reference=$1
probes=$2
p=$probes.`basename $reference .fa`
p=`basename $p .fa`
echo "naive probes: `grep -c ">" $probes`" > $p.log
# fine-grained BLAT to see where designed probes actually match
blat $reference $probes $p.blat -tileSize=10 -fine
# filter out BLAT matches that are imperfect, gapped, or partial
filter-blat-matches.py -R $p.blat > $p.filtered.dat
#
blat-probes.R.sh $p.filtered.dat $probes $p
cat $p.log
