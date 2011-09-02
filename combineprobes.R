# work in progress?
combineprobes=function(mapping,lvls,data){
	out=data.frame()
	cols=colnames(data)
	cols=cols[2:length(cols)]
	print(cols)
	print(lvls)
	for (l in 1:length(lvls)) {
		lvl=lvls[l]
		print(lvl)
		out[l,1]=lvl
		probes=mapping$probe[which(mapping$transcriptId==lvl)]
		for (p in 1:length(probes)) {
			prb=probes[p]
			print(prb)
			for (c in 1:length(cols)) {
				col=cols[c]
				print(col)
				i=grep(prb,data$probe)
				print(i)
				if(length(i)==1){
					val=data[[col]][i]
					out[l,c+1] = out[l,c+1]+val
				}
			}
		}
	}
	return(out)
}
