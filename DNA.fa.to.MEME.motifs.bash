for f in `cat $1`; do
	# hardcoded name filter
	prefix=`echo $f | sed -e 's#2e1c.P.OT3.FL11.##' | sed -e 's#/design_mutations.dna##'`
	echo $prefix
	# hardcoded start and end (rosetta dna .fa files read through both strands)
#	DNAMotif.py --bg halo.bg.file --MEME -p -o $prefix --start 1 --end 17 --pseudocounts 1 $f
	DNAMotif.py --bg halo.bg.file --MEME -p -o $prefix --start 2 --end 16 --pseudocounts 1 $f
done
