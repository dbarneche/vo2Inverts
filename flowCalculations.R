rm(list=ls())
source('paths.R')
source('R/functions-analyses.R')

flow            <-  readFile('data/fieldOxygenFlow.csv')
flow$flowCmSeg  <-  flow$Flow..m.s.1. * 100
round(tapply(flow$flowCmSeg, flow$Location, mean, na.rm=TRUE), 1)
round(tapply(flow$flowCmSeg, flow$Location, sd, na.rm=TRUE), 1)
