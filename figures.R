library(LoLinR)
library(plyr)
library(extrafont)
library(fontcm)
loadfonts()

rm(list=ls())
source('paths.R')
source('R/functions-figures.R')
load('output/RDatafiles/analyses.RData')

toPdf(comparisons(), 'output/figures/comparisons.pdf', width=9.5, height=4.5)
embed_fonts('output/figures/comparisons.pdf')

toPdf(comparisonsCIs(), 'output/figures/comparisonsCIs.pdf', width=9.5, height=4.5)
embed_fonts('output/figures/comparisonsCIs.pdf')

toPdf(comparisonsErect(), 'output/figures/comparisonsErect.pdf', width=6, height=6)
embed_fonts('output/figures/comparisonsErect.pdf')

toPdf(comparisonsErectCIs(), 'output/figures/comparisonsErectCIs.pdf', width=6, height=6)
embed_fonts('output/figures/comparisonsErectCIs.pdf')

toPdf(michaelisMentenJAGS(), 'output/figures/michaelisMentenJAGS.pdf', width=7.5, height=10)
embed_fonts('output/figures/michaelisMentenJAGS.pdf')

toPdf(massScaling(), 'output/figures/massScaling.pdf', width=7, height=7)
embed_fonts('output/figures/massScaling.pdf')
