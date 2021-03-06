######################
# AUXILLIARY FUNCTIONS
######################
linearRescale <- function(x, rOut) {
  p <- (x - min(x)) / (max(x) - min(x))
  rOut[[1]] + p * (rOut[[2]] - rOut[[1]])
}

violinPlotAndCumSumHist  <-  function(data, meanNative, meanInvasive) {
	d         <-  density(data$X.AS, from=min(data$X.AS), to=max(data$X.AS))
	y1        <-  d$x
	x1        <-  d$y
	x1        <-  linearRescale(x1, c(0, 0.3))
	negX      <-  unique(data$LocationNum) - x1
	polygon(c(negX, rep(max(negX), 2)), c(y1, y1[length(y1)], y1[1]), col=transparentColor('grey50', 0.6), border='grey50')
	text(max(negX), 160, unique(data$Location), adj=c(0.5,0), font=3)

	minBrk    <-  (floor((min(data$X.AS) + 5)/10)*10)-5
	maxBrk    <-  (ceiling((max(data$X.AS) + 5)/10)*10)-5
	h         <-  hist(data$X.AS, breaks=seq(minBrk,maxBrk,5), plot=FALSE)
	x2        <-  unique(data$LocationNum) + linearRescale(cumsum(h$counts), c(0,0.5))
	sapply(seq_along(x2), function(k, x2) {lines(c(x2[1], rep(x2[k], 2), rep(x2[1], 2)), c(rep(h$breaks[k], 2), rep(h$breaks[k]+5, 2), h$breaks[k]))}, x2=x2)
	percentageNative    <-  length(data$X.AS[data$X.AS <= meanNative])/nrow(data)
	percentageInvasive  <-  length(data$X.AS[data$X.AS <= meanInvasive])/nrow(data)
	cat(unique(data$LocationNum), 'native:', percentageNative, 'invasive:', percentageInvasive, '\n')
	lines(rep(unique(data$LocationNum)+percentageNative/2, 2), c(0, maxBrk), lty=2, col='dodgerblue')
	lines(rep(unique(data$LocationNum)+percentageInvasive/2, 2), c(0, maxBrk), lty=2, col='tomato')
}

#########
# FIGURES
#########
fig1  <-  function(o2tab, michaelisMentenFit) {
	par(mfrow=c(5,3), omi=c(1.1,1.1,0.1,0.1), mai=rep(0.3,4))
	d_ply(o2tab, .(species), function(x, fullModel, modelSummary) {
		species   <-  unique(x$sppNum)
		satRange  <-  seq(min(x$o2sat), max(x$o2sat), length.out = 50)
		asymp     <-  exp(modelSummary['lnA', 'mean'] + modelSummary[paste0('r[', species,',1]'), 'mean'])
		denPar    <-  exp(modelSummary['lnB', 'mean'] + modelSummary[paste0('r[', species,',2]'), 'mean'])
		plot(boundedO2Vol ~ o2sat, data = x, xlab = '', ylab = '', xlim = c(0,110), ylim = c(0,1.3), type = 'n', axes = FALSE)
		box()
		axis(1)
		axis(2, las = 1)
		points(boundedO2Vol ~ o2sat, data = x, pch = 16, col = transparentColor('dodgerblue2', 0.4))
		yline     <-  (asymp * satRange) / (denPar + satRange)
		vol100    <-  (asymp * 100) / (denPar + 100)
		sat50     <-  (vol100 * 0.5 * denPar) / (asymp - vol100 * 0.5)
		proportionalLabel(0.9, 0.9, substitute(dot('V')*'O'[2]*' @ 100%'==A, list(A = round(vol100, 2))), adj = c(1,0.5))
		lines(c(0, sat50), rep(vol100 * 0.5, 2), lty = 2, lwd = 1.5)
		lines(rep(sat50, 2), c(0, vol100 * 0.5), lty = 2, lwd = 1.5)
		arrows(sat50 + 23, y0 = vol100 * 0.25, x1 = sat50 + 3, length = 0.05, angle = 30, code = 2, lwd = 2)

		mcmcMat  <-  fullModel$BUGSoutput$sims.matrix
		for(i in 1:nrow(mcmcMat)) {
			mcmcAsymp     <-  exp(mcmcMat[i, 'lnA'] + mcmcMat[i, paste0('r[', species, ',1]')])
			mcmcDenPar    <-  exp(mcmcMat[i, 'lnB'] + mcmcMat[i, paste0('r[', species, ',2]')])
			mcmcYline     <-  (mcmcAsymp * satRange) / (mcmcDenPar + satRange)
			lines(satRange, mcmcYline, col = transparentColor('grey30', 0.05))
		}
		spName  <-  unique(x$species)
		genus   <-  strsplit(spName, ' ')[[1]][1]
		spec    <-  strsplit(spName, ' ')[[1]][2]
		lines(satRange, yline, col = 'tomato', lty = 1, lwd = 1.5)
		if(spec == 'sp')
			proportionalLabel(0.03, 1.1, substitute(italic(a) * ' ' * b, list(a = genus, b = spec)), adj = c(0,0.5), xpd = NA)
		else
			proportionalLabel(0.03, 1.1, substitute(italic(a) * ' ' * italic(b), list(a = genus, b = spec)), adj = c(0,0.5), xpd = NA)
	}, fullModel = michaelisMentenFit$mmfit, modelSummary = michaelisMentenFit$mmjagsout)
	mtext('Oxygen level (% air saturation)', side = 1, line = 1.5, outer = TRUE, cex = 1.3)
	mtext(substitute('Relative '*dot('V')*'O'[2]*' [0, 1]'), side = 2, line = 1, outer = TRUE, cex = 1.3)
}

fig2  <-  function(michaelisMentenFit, mcmcAnova) {
	mcmcMat  <-  michaelisMentenFit$mmfit$BUGSoutput$sims.matrix
	par(mfrow = c(1, 3), mai = c(1.02,0.8,0.82,0), omi = c(0,0.25,0,0.25), cex = 1, cex.axis = 1.2, cex.lab = 1.4, xpd = NA)
	x  <-  jitter(rep(c(1,2), each = nrow(mcmcMat)))
	y  <-  c(mcmcAnova$results$invasive, mcmcAnova$results$native)
	plot(x, y, axes = FALSE, type = 'n', xlab = 'Status', ylab = substitute('C'[CO[2]]*' (% air saturation)'), xlim = c(0.5,2.5), ylim = c(2, 20))
	box()
	axis(1, at = c(1,2), labels = c('Invasive', 'Native'))
	axis(2, at = seq(2, 20, 6), las = 1)
	points(x, y, pch = 16, col = transparentColor(rep(c('tomato', 'dodgerblue2'), each = nrow(mcmcMat)), 0.1), cex = 1.3)
	proportionalLabel(0.05, 0.95, '(a)', adj = c(0.5, 0.5), font = 2, cex = 1.1)

	par(mai = c(1.02,0.6,0.82,0.2))
	y  <-  c(mcmcAnova$results$erect, mcmcAnova$results$flat)
	plot(x, y, axes = FALSE, type = 'n', xlab = 'Shape', ylab = '', xlim = c(0.5,2.5), ylim = c(2, 20))
	box()
	axis(1, at = c(1,2), labels = c('Erect', 'Flat'))
	axis(2, las = 1, at = seq(2, 20, 6), labels = NA)
	points(x, y, pch = 16, col = transparentColor(rep(c('tomato', 'dodgerblue2'), each = nrow(mcmcMat)), 0.1), cex = 1.3)
	proportionalLabel(0.05, 0.95, '(b)', adj = c(0.5, 0.5), font = 2, cex = 1.1)

	x   <-  jitter(rep(c(1,2), each = nrow(mcmcMat)))
	cl  <-  transparentColor(rep(c('tomato', 'dodgerblue2'), each = nrow(mcmcMat)), 0.1)
	y   <-  c(mcmcAnova$resultsErect$invasive, mcmcAnova$resultsErect$native)
	par(mai = c(1.02,0.4,0.82,0.4))
	plot(x, y, axes = FALSE, type = 'n', xlab = 'Status', ylab = '', xlim = c(0.5,2.5), ylim = c(2, 20), xpd = NA)
	proportionalLabel(0.5, 1.05, 'Erect species', adj = c(0.5, 0.5), xpd = NA, font = 3, cex = 1.4)
	box()
	axis(1, at = c(1,2), labels = c('Invasive', 'Native'))
	axis(2, las = 1, at = seq(2, 20, 6), labels = NA)
	points(x, y, pch = 16, col = cl, cex = 1.3)
	proportionalLabel(0.05, 0.95, '(c)', adj = c(0.5, 0.5), font = 2, cex = 1.1)
}

fig3  <-  function(fieldFlow, o2tab, michaelisMentenFit) {
	## Calculate and plot the two histograms
	par(omi = rep(0.5, 4), cex = 1)
	plot(NA, xlab='', ylab='Oxygen level (% air sat.)', type='n', axes=FALSE, cex.lab=1.2, xpd=NA, xlim=c(0.5, 5.5), ylim=c(0,170), yaxs='i')
	box()
	axis(1, at=seq(1,5.5,0.5), labels=rep(c(0,1),5))
	axis(2, las=1)
	
	mcmcMat              <-  michaelisMentenFit$mmfit$BUGSoutput$sims.matrix
	status               <-  o2tab$status[match(1:14, o2tab$sppNum)]
	dat                  <-  data.frame()
	for(i in 1:nrow(mcmcMat)) {
		asymp       <-  exp(mcmcMat[i,'lnA'] + mcmcMat[i,paste0('r[', 1:14, ',1]')])
		denPar      <-  exp(mcmcMat[i,'lnB'] + mcmcMat[i,paste0('r[', 1:14, ',2]')])
		vol100      <-  (asymp * 100) / (denPar + 100)
		cpo2s5      <-  (vol100*0.95 * denPar) / (asymp - vol100*0.95)
		dat         <-  rbind(dat, data.frame(invasive = mean(cpo2s5[status == 'invasive']), native = mean(cpo2s5[status == 'native'])))
	}

	meanNative    <-  mean(dat$native)
	meanInvasive  <-  mean(dat$invasive)
	lapply(seq(1.5,4.5,1), function(x) {
		polygon(c(x+0.005, x+0.2, x+0.2, x+0.005, x+0.005), c(-5,-5, 5, 5, -5), col = 'white', border = NA, xpd = NA)
	})
	polygon(c(par('usr')[1]+0.005, 1-0.3, 1-0.3, par('usr')[1]+0.005, par('usr')[1]+0.005), c(-5,-5, 5, 5, -5), col = 'white', border = NA, xpd = NA)
	polygon(c(5.505, par('usr')[2]-0.005, par('usr')[2]-0.005, 5.505, 5.505), c(-5,-5, 5, 5, -5), col = 'white', border = NA, xpd = NA)

	d_ply(fieldFlow, .(LocationNum), violinPlotAndCumSumHist, meanNative, meanInvasive)
	lines(par('usr')[1:2], rep(meanNative, 2), lty = 2, col = 'dodgerblue2')
	lines(par('usr')[1:2], rep(meanInvasive, 2), lty = 2, col = 'tomato')
}
