library(R2jags)
library(plyr)

rm(list=ls())
load('output/RDatafiles/michaelisMentenFit.RData')

##########################
# GET VO2s AT CRITICAL PO2
# FOR EACH INDIVIDUAL
##########################
asymp   <-  exp(mmjagsout['lnA', 'mean'] + mmjagsout[paste0('r[', 1:14, ',1]'), 'mean'])
denPar  <-  exp(mmjagsout['lnB', 'mean'] + mmjagsout[paste0('r[', 1:14, ',2]'), 'mean'])
vol100  <-  (asymp * 100) / (denPar + 100)
cpo2s   <-  (vol100 * 0.5 * denPar) / (asymp - vol100 * 0.5) # same as cpo2s   <-  (-1 * denPar) / (1 - asymp / (vol100 * 0.5))
vo2AtCpo2  <-  ddply(o2tab, .(sppNum, column), function(x, asymp, denPar, cpo2s) {
	spid       <-  unique(x$sppNum)
	asymp      <-  asymp[spid]
	denPar     <-  denPar[spid]
	cpo2       <-  cpo2s[spid]
	vo2AtCpo2  <-  (asymp * max(x$o2volume) * cpo2) / (denPar + cpo2)
	data.frame(dry_mass_g = unique(x$dry_mass_g), vo2AtCpo2 = vo2AtCpo2, stringsAsFactors = FALSE)
}, asymp = asymp, denPar = denPar, cpo2s = cpo2s)

#########################################
# LOG-LOG MASS DEPENDENCE METABOLIC RATES
#########################################
set.seed(1)
mrJags       <-  list('lnO2volume' = log(vo2AtCpo2$vo2AtCpo2), 'id' = vo2AtCpo2$sppNum, 'lnDry_mass_g' = log(vo2AtCpo2$dry_mass_g))
mrinits1     <-  list(A = 0.5,  lnBo = -2, tauMu = 1)
mrinits2     <-  list(A = 0,    lnBo = 0,  tauMu = 3)
mrinits3     <-  list(A = -0.2, lnBo = -2, tauMu = 5)
mrfit        <-  jags(data=mrJags, inits=list(mrinits1, mrinits2, mrinits3), parameters.to.save=c('A', 'lnBo', 'varMu', 'varR', 'r'), model.file='massRates.bug', n.chains=3, n.iter=1.5e6, DIC=TRUE, n.thin=375)
mrfit        <-  autojags(mrfit, n.iter=1.5e6, n.thin=750, n.update=100)
mrjagsout    <-  mrfit$BUGSoutput$summary

rm(list=ls()[!(ls() %in% c('vo2AtCpo2', 'mrfit', 'mrjagsout'))])
save.image('output/RDatafiles/massRates.RData')
