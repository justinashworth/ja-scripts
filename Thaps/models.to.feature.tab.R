d=read.delim('tps.all.models.ja.tsv')
dd = data.frame( id=d$gid, type='CDS', name=d$tag, contig=d$genome, start_pos=d$start, end_pos=d$end, strand=d$strand, description=paste('[',d$name,']',' ',d$protein,sep=''), chrom_position=paste(d$start,'..',d$end), organism='Thalassiosira psuedonana',GeneID=d$gid)
tsv(dd,'feature.tab.new')
