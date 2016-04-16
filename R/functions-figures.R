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
	usr       <-  par('usr')
	rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
	whiteGrid()
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
		plot(boundedMassSpecificO2Vol ~ o2sat, data=x, xlab='', ylab='', xlim=c(0,110), ylim=c(0,1.3), type='n', axes=FALSE)
		usr       <-  par('usr')
		rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
		whiteGrid()
		box()
		axis(1)
		axis(2, las=1)
		points(boundedMassSpecificO2Vol ~ o2sat, data=x, pch=16, col=transparentColor('dodgerblue2', 0.4))
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
	mtext(substitute('Relative mass-specific '*dot('V')*'O'[2]*' [0, 1]'), side=2, line=1, outer=TRUE, cex=1.3)
}

comparisons  <-  function(mcmcMat = mmfit$BUGSoutput$sims.matrix) {
	par(mfrow=c(1, 2), mai=c(1.02,1.12,0.82,0.30), omi=c(0,0.25,0,0.25), cex.axis=1.2, cex.lab=1.4, xpd=NA)
	x  <-  jitter(rep(c(1,2), each=nrow(mcmcMat)))
	y  <-  c(results$invasive, results$native)
	plot(x, y, axes=FALSE, type='n', xlab='Status', ylab=substitute('Critical PO'[2]*' (%)'), xlim=c(0.5,2.5), ylim=c(5, 35))
	usr       <-  par('usr')
	rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
	whiteGrid()
	box()
	axis(1, at=c(1,2), labels=c('Invasive', 'Native'))
	axis(2, las=1)
	points(x, y, pch=16, col=transparentColor(rep(c('tomato', 'dodgerblue2'), each=nrow(mcmcMat)), 0.1), cex=1.3)

	y  <-  c(results$erect, results$flat)
	par(mai=c(1.02,0.52,0.82,0.90))
	plot(x, y, axes=FALSE, type='n', xlab='Shape', ylab='', xlim=c(0.5,2.5), ylim=c(5, 35))
	usr       <-  par('usr')
	rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
	whiteGrid()
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
	plot(x, y, axes=FALSE, type='n', xlab='Status', ylab=substitute('Critical PO'[2]*' (%)'%+-%' S.E.'), xlim=c(0.5,2.5), ylim=c(-5, 50))
	usr       <-  par('usr')
	rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
	whiteGrid()
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
	usr       <-  par('usr')
	rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
	whiteGrid()
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
	plot(x, y, axes=FALSE, type='n', xlab='Status', ylab=substitute('Critical PO'[2]*' (%)'%+-%' S.E.'), xlim=c(0.5,2.5), ylim=c(10, 50), xpd=NA)
	proportionalLabel(0.5, 1.1, 'Erect species', adj=c(0.5, 0.5), xpd=NA, font=3, cex=1.4)
	usr       <-  par('usr')
	rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
	whiteGrid()
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
	plot(x, y, axes=FALSE, type='n', xlab='Status', ylab=substitute('Critical PO'[2]*' (%)'%+-%' S.E.'), xlim=c(0.5,2.5), ylim=c(-5, 60), xpd=NA)
	proportionalLabel(0.5, 1.1, 'Erect species', adj=c(0.5, 0.5), xpd=NA, font=3, cex=1.4)
	usr       <-  par('usr')
	rect(usr[1], usr[3], usr[2], usr[4], col='grey90', border=NA)
	whiteGrid()
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
