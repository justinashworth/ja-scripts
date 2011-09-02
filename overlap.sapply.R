sapply(old.trh3.chr.lt.0.05[,1],function(x){TRUE %in% sapply(trh3.chr.lt.0.05[,1],function(y){abs(x-y)<1000})})
