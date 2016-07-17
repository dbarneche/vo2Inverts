######################
# AUXILLIARY FUNCTIONS
######################
toDev <- function(expr, dev, filename, ..., verbose=TRUE) {
  if ( verbose )
    cat(sprintf('Creating %s\n', filename))
  dev(filename, family='CM Roman', ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

toPdf <- function(expr, filename, ...) {
  toDev(expr, pdf, filename, ...)
}

rounded  <-  function(x, rounding=2) {
  format(round(x, rounding), nsmall=rounding)
}

linearRescale <- function(x, rOut) {
  p <- (x - min(x)) / (max(x) - min(x))
  rOut[[1]] + p * (rOut[[2]] - rOut[[1]])
}

#########
# FIGURES
#########
massScaling  <-  function() {
	spid   <-  read.csv('data/speciesID.csv', header=TRUE, stringsAsFactors=FALSE)
	maxO2  <-  daply(o2tab, .(column), function(x)mean(x$o2volume[x$o2sat >= quantile(x$o2sat, probs=0.99, type=2)]))
	par(omi=rep(0.5, 4), cex=1)
	x  <-  log(spid$dry_mass_g)
	y  <-  log(maxO2)
	plot(x, y, xlab='ln dry mass (mg)', ylab=substitute('ln '*dot('V')*'O'[2]*' (ml h'^{-1}*')'), type='n', axes=FALSE, cex.lab=1.2, xpd=NA)
	box()
	axis(1)
	axis(2, las=1)
	points(y ~ x, pch=16, col=transparentColor('dodgerblue2', 0.8), cex=1.5)
}

michaelisMentenJAGS  <-  function() {
	par(mfrow=c(5,3), omi=c(1.1,1.1,0.1,0.1), mai=rep(0.3,4))
	d_ply(o2tab, .(species), function(x, fullModel, modelSummary) {
		species   <-  unique(x$sppNum)
		satRange  <-  seq(min(x$o2sat), max(x$o2sat), length.out=50)
		asymp     <-  exp(modelSummary['lnA', 'mean'] + modelSummary[paste0('r[', species,',1]'), 'mean'])
		denPar    <-  exp(modelSummary['lnB', 'mean'] + modelSummary[paste0('r[', species,',2]'), 'mean'])
		plot(boundedO2Vol ~ o2sat, data=x, xlab='', ylab='', xlim=c(0,110), ylim=c(0,1.3), type='n', axes=FALSE)
		box()
		axis(1)
		axis(2, las=1)
		points(boundedO2Vol ~ o2sat, data=x, pch=16, col=transparentColor('dodgerblue2', 0.4))
		yline     <-  (asymp * satRange) / (denPar + satRange)
		vol100    <-  (asymp * 100) / (denPar + 100)
		sat50     <-  (vol100*0.5 * denPar) / (asymp - vol100*0.5)
		proportionalLabel(0.9, 0.9, substitute(dot('V')*'O'[2]*' @ 100%'==A, list(A=round(vol100, 2))), adj=c(1,0.5))
		lines(c(0, sat50), rep(vol100*0.5, 2), lty=2)
		lines(rep(sat50, 2), c(0, vol100*0.5), lty=2)
		mcmcMat  <-  fullModel$BUGSoutput$sims.matrix
		for(i in 1:nrow(mcmcMat)) {
			mcmcAsymp     <-  exp(mcmcMat[i, 'lnA'] + mcmcMat[i, paste0('r[', species, ',1]')])
			mcmcDenPar    <-  exp(mcmcMat[i, 'lnB'] + mcmcMat[i, paste0('r[', species, ',2]')])
			mcmcYline     <-  (mcmcAsymp * satRange) / (mcmcDenPar + satRange)
			lines(satRange, mcmcYline, col=transparentColor('grey30', 0.05))
		}
		spName  <-  unique(x$species)
		genus   <-  strsplit(spName, ' ')[[1]][1]
		spec    <-  strsplit(spName, ' ')[[1]][2]
		lines(satRange, yline, col='tomato', lty=1, lwd=1.5)
		if(spec == 'sp')
			proportionalLabel(0.03, 1.1, substitute(italic(a) * ' ' * b, list(a=genus, b=spec)), adj=c(0,0.5), xpd=NA)
		else
			proportionalLabel(0.03, 1.1, substitute(italic(a) * ' ' * italic(b), list(a=genus, b=spec)), adj=c(0,0.5), xpd=NA)
	}, fullModel=mmfit, modelSummary=mmjagsout)
	mtext(substitute('PO'[2]*' (%)'), side=1, line=1.5, outer=TRUE, cex=1.3)
	mtext(substitute('Relative '*dot('V')*'O'[2]*' [0, 1]'), side=2, line=1, outer=TRUE, cex=1.3)
}

comparisons  <-  function(mcmcMat = mmfit$BUGSoutput$sims.matrix) {
	par(mfrow=c(1, 2), mai=c(1.02,1.12,0.82,0.30), omi=c(0,0.25,0,0.25), cex.axis=1.2, cex.lab=1.4, xpd=NA)
	x  <-  jitter(rep(c(1,2), each=nrow(mcmcMat)))
	y  <-  c(results$invasive, results$native)
	plot(x, y, axes=FALSE, type='n', xlab='Status', ylab=substitute('Critical PO'[2]*' (%)'), xlim=c(0.5,2.5), ylim=c(5, 20))
	box()
	axis(1, at=c(1,2), labels=c('Invasive', 'Native'))
	axis(2, las=1)
	points(x, y, pch=16, col=transparentColor(rep(c('tomato', 'dodgerblue2'), each=nrow(mcmcMat)), 0.1), cex=1.3)

	y  <-  c(results$erect, results$flat)
	par(mai=c(1.02,0.52,0.82,0.90))
	plot(x, y, axes=FALSE, type='n', xlab='Shape', ylab='', xlim=c(0.5,2.5), ylim=c(5, 35))
	box()
	axis(1, at=c(1,2), labels=c('Erect', 'Flat'))
	axis(2, las=1, labels=NA)
	points(x, y, pch=16, col=transparentColor(rep(c('tomato', 'dodgerblue2'), each=nrow(mcmcMat)), 0.1), cex=1.3)
}

comparisonsCIs  <-  function(mcmcMat = mmfit$BUGSoutput$sims.matrix) {
	par(mfrow=c(1, 2), mai=c(1.02,1.12,0.82,0.30), omi=c(0,0.25,0,0.25), cex.axis=1.2, cex.lab=1.4, xpd=NA)
	x   <-  jitter(rep(c(1,2), each=nrow(mcmcMat)))
	cl  <-  transparentColor(rep(c('tomato', 'dodgerblue2'), each=nrow(mcmcMat)), 0.1)
	y   <-  c(results$invasive, results$native)
	y1  <-  y - c(results$invasive_se, results$native_se)
	y2  <-  y + c(results$invasive_se, results$native_se)
	plot(x, y, axes=FALSE, type='n', xlab='Status', ylab=substitute('Critical PO'[2]*' (%)'%+-%' S.E.'), xlim=c(0.5,2.5), ylim=c(5, 20))
	box()
	axis(1, at=c(1,2), labels=c('Invasive', 'Native'))
	axis(2, las=1)
	for(i in seq_along(x)) {
		lines(rep(x[i], 2), c(y1[i], y2[i]), pch=16, col=cl[i])
	}
	points(x, y, pch=21, col=transparentColor('grey30', 0.2), bg=cl, cex=1.3)
	av     <-  tapply(y, rep(c(1,2), each=nrow(mcmcMat)), median)
	avLow  <-  tapply(y1, rep(c(1,2), each=nrow(mcmcMat)), median)
	avHi   <-  tapply(y2, rep(c(1,2), each=nrow(mcmcMat)), median)
	lines(c(1,1), c(avLow[1], avHi[1])) 
	lines(c(2,2), c(avLow[2], avHi[2])) 
	points(c(1,2), av, pch=16)

	y   <-  c(results$erect, results$flat)
	y1  <-  y - c(results$erect_se, results$flat_se)
	y2  <-  y + c(results$erect_se, results$flat_se)
	par(mai=c(1.02,0.52,0.82,0.90))
	plot(x, y, axes=FALSE, type='n', xlab='Shape', ylab='', xlim=c(0.5,2.5), ylim=c(-5, 50))
	box()
	axis(1, at=c(1,2), labels=c('Erect', 'Flat'))
	axis(2, las=1, labels=NA)
	for(i in seq_along(x)) {
		lines(rep(x[i], 2), c(y1[i], y2[i]), pch=16, col=cl[i])
	}
	points(x, y, pch=21, col=transparentColor('grey30', 0.2), bg=cl, cex=1.3)
	av     <-  tapply(y, rep(c(1,2), each=nrow(mcmcMat)), median)
	avLow  <-  tapply(y1, rep(c(1,2), each=nrow(mcmcMat)), median)
	avHi   <-  tapply(y2, rep(c(1,2), each=nrow(mcmcMat)), median)
	lines(c(1,1), c(avLow[1], avHi[1])) 
	lines(c(2,2), c(avLow[2], avHi[2])) 
	points(c(1,2), av, pch=16)
}

comparisonsErect  <-  function(mcmcMat = mmfit$BUGSoutput$sims.matrix) {
	par(cex=1, omi=rep(0.5, 4), cex.axis=1.2, cex.lab=1.4)
	x   <-  jitter(rep(c(1,2), each=nrow(mcmcMat)))
	cl  <-  transparentColor(rep(c('tomato', 'dodgerblue2'), each=nrow(mcmcMat)), 0.1)
	y   <-  c(resultsErect$invasive, resultsErect$native)
	plot(x, y, axes=FALSE, type='n', xlab='Status', ylab=substitute('Critical PO'[2]*' (%)'%+-%' S.E.'), xlim=c(0.5,2.5), ylim=c(5, 20), xpd=NA)
	proportionalLabel(0.5, 1.1, 'Erect species', adj=c(0.5, 0.5), xpd=NA, font=3, cex=1.4)
	box()
	axis(1, at=c(1,2), labels=c('Invasive', 'Native'))
	axis(2, las=1)
	points(x, y, pch=16, col=cl, cex=1.3)
}

comparisonsErectCIs  <-  function(mcmcMat = mmfit$BUGSoutput$sims.matrix) {
	par(cex=1, omi=rep(0.5, 4), cex.axis=1.2, cex.lab=1.4)
	x   <-  jitter(rep(c(1,2), each=nrow(mcmcMat)))
	cl  <-  transparentColor(rep(c('tomato', 'dodgerblue2'), each=nrow(mcmcMat)), 0.1)
	y   <-  c(resultsErect$invasive, resultsErect$native)
	y1  <-  y - c(resultsErect$invasive_se, resultsErect$native_se)
	y2  <-  y + c(resultsErect$invasive_se, resultsErect$native_se)
	plot(x, y, axes=FALSE, type='n', xlab='Status', ylab=substitute('Critical PO'[2]*' (%)'%+-%' S.E.'), xlim=c(0.5,2.5), ylim=c(5, 20), xpd=NA)
	proportionalLabel(0.5, 1.1, 'Erect species', adj=c(0.5, 0.5), xpd=NA, font=3, cex=1.4)
	box()
	axis(1, at=c(1,2), labels=c('Invasive', 'Native'))
	axis(2, las=1)
	for(i in seq_along(x)) {
		lines(rep(x[i], 2), c(y1[i], y2[i]), pch=16, col=cl[i])
	}
	points(x, y, pch=21, col=transparentColor('grey30', 0.2), bg=cl, cex=1.3)
	av     <-  tapply(y, rep(c(1,2), each=nrow(mcmcMat)), median)
	avLow  <-  tapply(y1, rep(c(1,2), each=nrow(mcmcMat)), median)
	avHi   <-  tapply(y2, rep(c(1,2), each=nrow(mcmcMat)), median)
	lines(c(1,1), c(avLow[1], avHi[1])) 
	lines(c(2,2), c(avLow[2], avHi[2])) 
	points(c(1,2), av, pch=16)
}

plotViolin  <-  function(data, relativeDensCeiling = 0.3) {
	dens    <-  density(data$X.AS, from=min(data$X.AS), to=max(data$X.AS))
	quants  <-  quantile(dens$x, probs = c(0.025, 0.975), type = 2)
	y       <-  dens$y[dens$x >= quants[1] & dens$x <= quants[2]]
	x       <-  dens$x[dens$x >= quants[1] & dens$x <= quants[2]]
	y       <-  linearRescale(y, c(0, relativeDensCeiling))
	negY    <-  unique(data$LocationNum) - y
	posY    <-  unique(data$LocationNum) + y
	polygon(c(negY, posY[length(posY):1], negY[1]), c(x, x[length(x):1], x[1]), col=transparentColor('dodgerblue2', 0.6), border=NA)
	lines(c(negY, posY[length(posY):1], negY[1]), c(x, x[length(x):1], x[1]), col='dodgerblue2')
	text(unique(data$LocationNum), max(data$X.AS), substitute(a %+-% b, list(a=rounded(mean(data$Flow..m.s.1., na.rm=TRUE), 2), b=rounded(sd(data$Flow..m.s.1., na.rm=TRUE), 2))), adj=c(0.5, 0), cex = 0.8, font = 3)
}

violinsPlot  <-  function() {
	par(omi = rep(0.5, 4), cex = 1)
	plot(NA, xlab='', ylab='Oxygen level (% air sat.)', type='n', axes=FALSE, cex.lab=1.2, xpd=NA, xlim=c(0.5, 5.5), ylim=range(fieldFlow$X.AS))
	box()
	axis(1, at=c(1:5), labels=sort(unique(fieldFlow$Location)), las=3)
	axis(2, las=1)
	d_ply(fieldFlow, .(Location), plotViolin)
}


fig1  <-  function() {
	## Calculate and plot the two histograms
	par(omi = rep(0.5, 4), cex = 1)
	plot(NA, xlab='', ylab='Oxygen level (% air sat.)', type='n', axes=FALSE, cex.lab=1.2, xpd=NA, xlim=c(0.5, 5.5), ylim=c(0,160), yaxs='i')
	box()
	axis(1, at=seq(1,5.5,0.5), labels=rep(c(0,1),5))
	axis(2, las=1)
	meanNative    <-  mean(results$native)
	meanInvasive  <-  mean(results$invasive)
	lapply(seq(1.5,4.5,1), function(x) {
		polygon(c(x+0.005, x+0.2, x+0.2, x+0.005, x+0.005), c(-5,-5, 5, 5, -5), col='white', border=NA, xpd=NA)
	})
	polygon(c(par('usr')[1]+0.005, 1-0.3, 1-0.3, par('usr')[1]+0.005, par('usr')[1]+0.005), c(-5,-5, 5, 5, -5), col='white', border=NA, xpd=NA)
	polygon(c(5.505, par('usr')[2]-0.005, par('usr')[2]-0.005, 5.505, 5.505), c(-5,-5, 5, 5, -5), col='white', border=NA, xpd=NA)

	d_ply(fieldFlow, .(LocationNum), violinPlotAndCumSumHist, meanNative, meanInvasive)
	lines(par('usr')[1:2], rep(meanNative, 2), lty=2, col='dodgerblue2')
	lines(par('usr')[1:2], rep(meanInvasive, 2), lty=2, col='tomato')
}

violinPlotAndCumSumHist  <-  function(data, meanNative, meanInvasive) {
	d         <-  density(data$X.AS, from=min(data$X.AS), to=max(data$X.AS))
	y1        <-  d$x
	x1        <-  d$y
	x1        <-  linearRescale(x1, c(0, 0.3))
	negX      <-  unique(data$LocationNum) - x1
	polygon(c(negX, rep(max(negX), 2)), c(y1, y1[length(y1)], y1[1]), col=transparentColor('grey50', 0.6), border='grey50')
	text(max(negX), y1[length(y1)]+10, unique(data$Location), adj=c(0.5,0), font=3)

	minBrk    <-  (floor((min(data$X.AS) + 5)/10)*10)-5
	maxBrk    <-  (ceiling((max(data$X.AS) + 5)/10)*10)-5
	h         <-  hist(data$X.AS, breaks=seq(minBrk,maxBrk,5), plot=FALSE)
	x2        <-  unique(data$LocationNum) + linearRescale(cumsum(h$counts), c(0,0.5))
	for(k in seq_along(x2)) {
		lines(c(x2[1], rep(x2[k], 2), rep(x2[1], 2)), c(rep(h$breaks[k], 2), rep(h$breaks[k]+5, 2), h$breaks[k]))
	}
	percentageNative    <-  length(data$X.AS[data$X.AS <= meanNative])/nrow(data)
	percentageInvasive  <-  length(data$X.AS[data$X.AS <= meanInvasive])/nrow(data)
	lines(rep(unique(data$LocationNum)+percentageNative, 2), c(0, maxBrk), lty=2, col='dodgerblue')
	lines(rep(unique(data$LocationNum)+percentageInvasive, 2), c(0, maxBrk), lty=2, col='tomato')
}

