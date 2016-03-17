library(plyr)
set.seed(1)

rm(list=ls())
# load volume data and air saturation data
vol  <-  read.csv('data/vol_2015_12_17.csv', header=TRUE, stringsAsFactors=FALSE)
ast  <-  read.csv('data/airSat_2015_12_17.csv', header=TRUE, stringsAsFactors=FALSE)

# make sure that the columns match in size and structure
dim(vol) == dim(ast)

# make sure that NAs are placed correctly
for(i in 1:ncol(vol)) {
  cat(i, ': ', length(vol[!is.na(vol[,i]),i]) == length(ast[!is.na(ast[,i]),i]), '\n')
}; rm(i)
    
# check quality of data
#dim(data)
#all(sapply(data, class) == 'numeric')
#head(data)

# volume data - standardise columns by their respective maxima, i.e. [0,1]
for(i in 1:ncol(vol)) {
  vol[, i]  <-  vol[,i]/max(vol[,i], na.rm=TRUE)
  cat(range(vol[,i], na.rm=TRUE), '\n')
}; rm(i)

# create functions
source('R/functions.R')
pco2s  <-  vector(mode='list', length=ncol(vol))
for(i in 1:ncol(vol)) {
  pco2s[[i]]  <-  pco2(po2=ast[,i], vo2=vol[,i], index=i)
}
pco2s   <-  do.call(rbind.data.frame, pco2s)
best    <-  ddply(pco2s, .(column), function(x) {
	x  <-  x[complete.cases(x), ]
	x[which.min(x$AICs),]
})
write.csv(pco2s, 'output/data/completeAICs.csv', row.names=FALSE)
write.csv(best, 'output/data/bestAICs.csv', row.names=FALSE)
