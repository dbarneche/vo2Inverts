library(plyr)
rm(list=ls())
source('database.R')

flow  <-  readFile('data/fieldOxygenFlow.csv')

table1  <-  ddply(o2tab, .(species), function(x) {
	data.frame('Growth shape' = unique(x$shape),
			   'Status' = unique(x$status),
			   'n' = length(unique(x$column)), stringsAsFactors=FALSE)
})

table2  <-  ddply(flow, .(LocationNum), function(x) {
	data.frame(Site=unique(x$Location),
		Mean=round(mean(x$X.AS, na.rm=TRUE), 2),
		SD=round(sd(x$X.AS, na.rm=TRUE), 2),
		Min=round(min(x$X.AS, na.rm=TRUE), 2),
		Max=round(max(x$X.AS, na.rm=TRUE), 2), stringsAsFactors=FALSE)
})
