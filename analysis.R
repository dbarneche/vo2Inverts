library(plyr)

rm(list=ls())
load('output/RDatafiles/michaelisMentenFit.RData')

##############################
# EXTRACT PARAMETERS FROM EACH
# MCMC SAMPLE AND FIT ANOVAS
# ALSO GET STANDARD ERRORS
##############################
mcmcMat              <-  mmfit$BUGSoutput$sims.matrix
species              <-  o2tab$species[match(1:14, o2tab$sppNum)]
shapes               <-  o2tab$shape[match(1:14, o2tab$sppNum)]
status               <-  o2tab$status[match(1:14, o2tab$sppNum)]
results              <-  data.frame(matrix(0, nrow(mcmcMat), 8))
names(results)       <-  c('invasive', 'invasive_se', 'native', 'native_se', 'flat', 'flat_se', 'erect', 'erect_se')
resultsErect         <-  data.frame(matrix(0, nrow(mcmcMat), 4))
names(resultsErect)  <-  c('invasive', 'invasive_se', 'native', 'native_se')

for(i in 1:nrow(mcmcMat)) {
	asymp       <-  exp(mcmcMat[i,'lnA'] + mcmcMat[i,paste0('r[', 1:14, ',1]')])
	denPar      <-  exp(mcmcMat[i,'lnB'] + mcmcMat[i,paste0('r[', 1:14, ',2]')])
	vol100      <-  (asymp * 100) / (denPar + 100)
	cpo2s       <-  (vol100*0.5 * denPar) / (asymp - vol100*0.5)
	cpo2sErect  <-  cpo2s[shapes == 'erect']
	statErect   <-  status[shapes == 'erect']
	mod1        <-  lm(cpo2s ~ status - 1)
	mod2        <-  lm(cpo2s ~ shapes - 1)
	mod3        <-  lm(cpo2sErect ~ statErect - 1)
	
	results$invasive[i]          <-  coef(mod1)['statusinvasive']
	results$invasive_se[i]       <-  coef(summary(mod1))['statusinvasive', 'Std. Error']
	results$native[i]            <-  coef(mod1)['statusnative']
	results$native_se[i]         <-  coef(summary(mod1))['statusnative', 'Std. Error']
	results$erect[i]             <-  coef(mod2)['shapeserect']
	results$erect_se[i]          <-  coef(summary(mod2))['shapeserect', 'Std. Error']
	results$flat[i]              <-  coef(mod2)['shapesflat']
	results$flat_se[i]           <-  coef(summary(mod2))['shapesflat', 'Std. Error']
	resultsErect$invasive[i]     <-  coef(mod3)['statErectinvasive']
	resultsErect$invasive_se[i]  <-  coef(summary(mod3))['statErectinvasive', 'Std. Error']
	resultsErect$native[i]       <-  coef(mod3)['statErectnative']
	resultsErect$native_se[i]    <-  coef(summary(mod3))['statErectnative', 'Std. Error']	
}

##############################
# SUMMARY AVERAGES FOR FLAT
# SPECIES
##############################
flatNative           <-  shapes == 'flat' & status == 'native'
flatInvasive         <-  shapes == 'flat' & status == 'invasive'
cpo2sFlat            <-  adply(mcmcMat, 1, function(x, flatNative, flatInvasive) {
	asymp       <-  exp(x['lnA'] + x[paste0('r[', 1:14, ',1]')])
	denPar      <-  exp(x['lnB'] + x[paste0('r[', 1:14, ',2]')])
	vol100      <-  (asymp * 100) / (denPar + 100)
	allPco2s    <-  (vol100*0.5 * denPar) / (asymp - vol100*0.5)
	data.frame(meanNative=mean(allPco2s[flatNative]), sdNative=sd(allPco2s[flatNative]), meanInvasive=mean(allPco2s[flatInvasive]), sdInvasive=sd(allPco2s[flatInvasive]))
}, flatNative=flatNative, flatInvasive=flatInvasive, .id=NULL)
cpo2sFlat            <-  apply(cpo2sFlat, 2, median)
	
rm(list=ls()[!(ls() %in% c('o2tab', 'mmfit', 'mmjagsout', 'results', 'resultsErect', 'cpo2sFlat'))])
save.image('output/RDatafiles/analyses_resultsAirSat@v50.RData')
