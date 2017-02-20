readFile  <-  function(filepath){
    read.csv(filepath, header=TRUE, stringsAsFactors=FALSE, na.strings=c('','NA'), strip.white=TRUE)  
}

cleanData  <-  function(vol, ast, spid) {
    o2tab  <-  data.frame()
    for(i in seq_len(ncol(vol))) {
        dat  <-  data.frame(species     =  spid$species[spid$column == i], 
                            dry_mass_g  =  spid$dry_mass_g[spid$column == i], 
                            status      =  tolower(spid$status[spid$column == i]), 
                            shape       =  tolower(spid$shape[spid$column == i]), 
                            column      =  i, 
                            o2volume    =  vol[!is.na(vol[,i]),i], 
                            o2sat       =  ast[!is.na(ast[,i]),i],
                            stringsAsFactors = FALSE)
        
        dat$boundedO2Vol  <-  dat$o2volume / max(dat$o2volume)
        o2tab             <-  rbind(o2tab, dat)
    }
    o2tab$sppNum  <-  as.numeric(as.factor(tolower(gsub(' ', '', o2tab$species))))
    o2tab
}

fitMichaelisMentenInJAGS  <-  function(o2tab, modelFile) {
    set.seed(1)
    mmJags       <-  list('boundedO2Vol' = o2tab$boundedO2Vol, 'id' = o2tab$sppNum, 'o2sat' = o2tab$o2sat)
    mminits1     <-  list(lnA = 0.5,  lnB = 2, tauMu = 1)
    mminits2     <-  list(lnA = 0,    lnB = 5, tauMu = 3)
    mminits3     <-  list(lnA = -0.2, lnB = 3, tauMu = 5)
    mmfit        <-  jags(data = mmJags, inits = list(mminits1, mminits2, mminits3), parameters.to.save = c('lnA', 'lnB', 'varMu', 'varR', 'r'), model.file = modelFile, n.chains = 3, n.iter = 1.5e6, DIC = TRUE, n.thin = 375)
    mmfit        <-  autojags(mmfit, n.iter = 1.5e6, n.thin = 750, n.update = 100)
    list(mmfit = mmfit, mmjagsout = mmfit$BUGSoutput$summary)
}

extractMcmcAnova  <-  function(o2tab, michaelisMentenFit) {
    ##############################
    # EXTRACT PARAMETERS FROM EACH
    # MCMC SAMPLE AND FIT ANOVAS
    # ALSO GET STANDARD ERRORS
    ##############################
    mcmcMat              <-  michaelisMentenFit$mmfit$BUGSoutput$sims.matrix
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
        cpo2s       <-  (vol100 * 0.5 * denPar) / (asymp - vol100 * 0.5) # same as cpo2s   <-  (-1 * denPar) / (1 - asymp / (vol100 * 0.5))
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
    list(results = results, resultsErect = resultsErect)
}

getAvCpo2sFlat  <-  function(michaelisMentenFit) {
    ##############################
    # SUMMARY AVERAGES FOR FLAT
    # SPECIES
    ##############################
    mcmcMat              <-  michaelisMentenFit$mmfit$BUGSoutput$sims.matrix
    flatNative           <-  shapes == 'flat' & status == 'native'
    flatInvasive         <-  shapes == 'flat' & status == 'invasive'
    cpo2sFlat            <-  adply(mcmcMat, 1, function(x, flatNative, flatInvasive) {
        asymp       <-  exp(x['lnA'] + x[paste0('r[', 1:14, ',1]')])
        denPar      <-  exp(x['lnB'] + x[paste0('r[', 1:14, ',2]')])
        vol100      <-  (asymp * 100) / (denPar + 100)
        allPco2s    <-  (vol100*0.5 * denPar) / (asymp - vol100*0.5)
        data.frame(meanNative    =  mean(allPco2s[flatNative]),
                   sdNative      =  sd(allPco2s[flatNative]), 
                   meanInvasive  =  mean(allPco2s[flatInvasive]), 
                   sdInvasive    =  sd(allPco2s[flatInvasive]))
    }, flatNative = flatNative, flatInvasive = flatInvasive, .id = NULL)
    apply(cpo2sFlat, 2, median)
}
