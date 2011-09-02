# "d" is tab-delimited BLAT table
# this adds
blattrack = function(d,name="track",ds=getDatasetDescription(),style="qseg") {
	track = data.frame(sequence=d$Tname,strand=d$strand,start=d$Tstart,end=d$Tend,value=1)
	addTrack(ds,track,name)
}
