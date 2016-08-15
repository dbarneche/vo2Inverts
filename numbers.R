rm(list=ls())
source('R/functions-analyses.R')
load('output/RDatafiles/analyses_resultsAirSat@v50.RData')

flatErectFoldChange  <-  1 / (results$flat / results$erect)
mean(flatErectFoldChange)
quantile(flatErectFoldChange, probs=c(0.025, 0.975), type=2)

invasiveNativeFoldChange  <-  1 / (results$invasive / results$native)
mean(invasiveNativeFoldChange)
quantile(invasiveNativeFoldChange, probs=c(0.025, 0.975), type=2)
