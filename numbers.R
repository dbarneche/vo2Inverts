library(plyr)

rm(list=ls())
source('R/functions-analyses.R')
load('output/RDatafiles/analyses_resultsAirSat@v50.RData')

flatErectFoldChange  <-  1 / (results$flat / results$erect)
mean(flatErectFoldChange)
quantile(flatErectFoldChange, probs=c(0.025, 0.975), type=2)

invasiveNativeFoldChange  <-  1 / (results$invasive / results$native)
mean(invasiveNativeFoldChange)
quantile(invasiveNativeFoldChange, probs=c(0.025, 0.975), type=2)

load('output/RDatafiles/massRates.RData')
mcmcMatBnots  <-  mrfit$BUGSoutput$sims.matrix
status        <-  o2tab$status[match(1:14, o2tab$sppNum)]

vo2sAtCpo2sAtOneGram  <-  t(apply(mcmcMatBnots, 1, function(x) {
	exp(x['lnBo'] + x[paste0('r[', 1:14, ',1]')])
}))

sppVO2Table  <-  adply(vo2sAtCpo2sAtOneGram, 2, function(x) {
	cis95  <-  quantile(x, probs = c(.025, .975), type = 2)
	data.frame(mean = mean(x), lower = cis95[1], upper = cis95[2], stringsAsFactors = FALSE)
}, .id = 'sppNum')

statusVO2McMcTable  <-  as.matrix(adply(vo2sAtCpo2sAtOneGram, 1, function(x, status) {
	means  <-  tapply(x, status, mean)
	data.frame(invasive = means[1], native = means[2], stringsAsFactors = FALSE)
}, .id = NULL, status = status))

statusVO2Table  <-  adply(statusVO2McMcTable, 2, function(x) {
	cis95  <-  quantile(x, probs = c(.025, .975), type = 2)
	data.frame(mean = mean(x), lower = cis95[1], upper = cis95[2], stringsAsFactors = FALSE)
}, .id = 'status')
