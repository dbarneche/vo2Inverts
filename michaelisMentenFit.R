library(R2jags)

rm(list=ls())
source('paths.R')
source('database.R')

#################################
# MICHAELIS-MENTEN MODEL IN JAGS
#################################
set.seed(1)
mmJags       <-  list('boundedO2Vol'=o2tab$boundedO2Vol, 'id'=o2tab$sppNum, 'o2sat'=o2tab$o2sat)
mminits1     <-  list(lnA = 0.5,  lnB = 2, tauMu=1)
mminits2     <-  list(lnA = 0,    lnB = 5, tauMu=3)
mminits3     <-  list(lnA = -0.2, lnB = 3, tauMu=5)
mmfit        <-  jags(data=mmJags, inits=list(mminits1, mminits2, mminits3), parameters.to.save=c('lnA', 'lnB', 'varMu', 'varR', 'r'), model.file='michaelisMenten.bug', n.chains=3, n.iter=1.5e6, DIC=TRUE, n.thin=375)
mmfit        <-  autojags(mmfit, n.iter=1.5e6, n.thin=750, n.update=100)
mmjagsout    <-  mmfit$BUGSoutput$summary #jags output

rm(list=ls()[!(ls() %in% c('o2tab', 'mmfit', 'mmjagsout'))])
save.image('output/RDatafiles/michaelisMentenFit.RData')
