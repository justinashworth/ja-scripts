# really basic script to combine two data matrices

# this is an existing data matrix, in tsv format with an 'id' column
# for example, an existing log2 ratios matrix
existing_data = 'CO2.trans.tsv'

# some new data matrix; needs an 'id' column that corresponds to the 'id' column in the existing data (it does not need to be sorted)
# data type should match existing data (e.g. log2 ratios to some internal control or other suitable control--best done separately, and before merging)
new_data = 'CO2.stead.tsv'

# R's merge() will nicely paste and reorganize the new and existing data together based on the 'id' column
# the 'all' flag means you get the union of the two datasets (default is the intersection, dropping any ids that do not occur in both datasets; all.x=T means don't drop any existing data)
merged_data = merge( read.delim(existing_data,as.is=T), read.delim(new_data,as.is=T), by='id', all=T)

# optional: remove any rows with all NA's
nnas = apply(merged_data,1,function(x){length(which(is.na(x[colnames(merged_data)!='id'])))})
merged_data = merged_data[nnas==0,]

# write a tsv
write.table(merged_data,'mergeddata.tsv',row.names=F,quote=F,sep='\t')

# repeat as necessary. (Alternatively, code up something more complicated that merges multiple files in one go.)

# below is optional stuff

# write an R table/matrix
mat = as.matrix(merged_data[,2:ncol(merged_data)])
write.table(mat,'mergeddata.tab')

# normalization: the appropriate normalization approach depends on the downstream application and statistical requirements (or lack thereof) for various statistical testing models

# example: of normalizing within-experiment ratios such that the standard deviations of all columns are the same:
mat[,] = apply(mat, 2, function(x){ x/sd(x,na.rm=T) })
write.table(mat,'mergeddata.colsdnorm.tab')

# example: further normalize such that the column means are all zero (effectively, data are now within-experiment Z-scores)
mat[,] = apply(mat, 2, function(x){ x - mean(x,na.rm=T) })
write.table(mat,'mergeddata.colZnorm.tab')

# other normalization approaches exist, e.g. quantile normalization, row normalization, MAD normalization...
