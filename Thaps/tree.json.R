# for writing a json tree files,
# e.g. for d3 'treemap' http://bl.ocks.org/mbostock/4063582#index.html
# steps:
# 1. select a dendrogram,
# 2. combine distal leaves (clusters):
#    cut() dendrogram at the correct height (as determined by dendextendRcpp_heights_per_k.dendrogram, e.g. height for k=500), keep $upper
# 3. break top of tree into subgroups: use cut() again on $upper, this time at an hpk corresponding to desired number of subtrees (e.g. height for k=20); keep $lower
# 4. run this script, write to file='tree.json'
tree_json = function(dnd,maxlevel=5,level=0,prefix='',file='',k="500"){
	tab=paste(rep(' ',level),collapse='')
	if(is.leaf(dnd)){
		cl = gsub('Branch ','',labels(dnd))
		if(exists('clids')) sz = length(clids[[k]][[as.integer(cl)]])
		else sz = 1
		tx=sprintf('%s{"name": "%s", "size": %i}',tab,cl,sz)
		cat(tx,append=T,file=file)
	}else if(level>=maxlevel){
		for(i in 1:length(labels(dnd))){
			cl = gsub('Branch ','',labels(dnd)[i])
			if(exists('clids')) sz = length(clids[[k]][[as.integer(cl)]])
			else sz = 1
			tx=sprintf('%s{"name": "%s", "size": %i}',tab,cl,sz)
			if(i!=length(labels(dnd))) tx = paste(tx,',',sep='')
			cat(tx,'\n',append=T,file=file)
		}
	}else{
#		if(length(labels(dnd))>0){
			tx=sprintf('%s{"name": "Branch%s", "children": [\n',tab,prefix)
			cat(tx,append=T,file=file)
			level = level + 1
			for(i in 1:length(dnd)){
				tree_json(dnd[[i]],maxlevel,file=file,level,paste(prefix,i,sep='_'))
				if(i!=length(dnd)) cat(',\n',append=T,file=file)
			}
			tx = sprintf('%s]}\n',tab)
			cat(tx,append=T,file=file)
#		}
	}
}

system('rm tree.json')
k = "400"
h = hpk[names(hpk)==k]
up = cut(as.dendrogram(hc),h=h)$upper
groups = "20"
h = hpk[names(hpk)==groups]
mid = cut(as.dendrogram(up),h=h)$lower
tree_json(mid,maxlevel=4,file='tree.json',k=k)
