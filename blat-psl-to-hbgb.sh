#!/bin/bash
outf=`basename $1 .psl`.hbgb
#tail -n +6 $1 | awk '{OFS="\t"}{print $10,$14,$9,$16,$17}' | sed "s/chr_\([0-9]*[a-z]\)_[0-9]*/chr_\1/" | sed "s/chr_/Chr/" > $outf
tail -n +6 $1 | awk '{OFS="\t"}{print $14,$9,$16,$17,$10}' | sed "s/chr_\([0-9]*[a-z]\)_[0-9]*/chr_\1/" | sed "s/chr_/Chr/" > $outf
#add-value-column.py -v $2 $outf
