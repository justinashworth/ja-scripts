#!/bin/bash
mast_exec=$HOME/bin/mast4.5
args=""
args="-bfile halo.bg.file -w -mt 0.001"
args="-bfile halo.bg.file"
#args="-bfile halo.bg.file -w"
for lrp in `cat lrps`; do
	motif_file="$lrp.MEME"
#	for tgt in *tgts.fa; do
#		$mast_exec $motif_file $tgt $args -hit_list -stdout > $lrp.promoter.hits
#		$mast_exec $motif_file $tgt
#	done
	# html output
	$mast_exec $motif_file all.tgt.seqs $args -oc mastout.$lrp
	# have to run again to get hit_list files (?)
	echo "target	strand	start	end	score	pval" > mast.$lrp
	$mast_exec $motif_file all.tgt.seqs $args -oc mastout.$lrp -hit_list | sed -E 's/ +/	/g' | grep -v "#" >> mast.$lrp
done
