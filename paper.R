library(knitr)
library(bibtex)
library(knitcitations)
library(plyr)

rm(list=ls())
source('R/functions-analyses.R')
load('output/RDatafiles/analyses_resultsAirSat@v50.RData')

# to .md
knit('text/lagos_et_al.Rmd', output=file.path(getwd(), 'text/lagos_et_al.md'), quiet=TRUE, encoding = 'utf-8')
knit('text/supplement.Rmd', output=file.path(getwd(), 'text/supplement.md'), quiet=TRUE, encoding = 'utf-8')

# to .docx
system('pandoc -o text/lagos_et_al.docx text/lagos_et_al.md -s -S --bibliography library.bib --csl ecology.csl')
system('pandoc -o text/supplement.docx text/supplement.md -s -S --bibliography library.bib --csl ecology.csl')
