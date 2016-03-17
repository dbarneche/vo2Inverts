michaelisMenten <- function(a, b, po2) 
  a * po2 / (b + po2)
michaelisMenten <- deriv(body(michaelisMenten), namevec = c("a", "b"), func = michaelisMenten)

powerFct <- function(a, b, po2) 
  a * po2 ^ b
powerFct <- deriv(body(powerFct), namevec = c("a", "b"), func = powerFct)

hyperbola <- function(a, b, c, po2) 
  a * po2 / (b + po2) + c
hyperbola <- deriv(body(hyperbola), namevec = c("a", "b", "c"), func = hyperbola)

pareto <- function(a, b, po2) 
  1 - (a * po2) ^ b
pareto <- deriv(body(pareto), namevec = c("a", "b"), func = pareto)

runAllModels  <-  function(po2, vo2) {
  dat  <-  data.frame(po2, vo2)
  dat  <-  dat[complete.cases(dat), ]
  dat[dat <= 0]  <-  1e-10
  micmen  <-  hyper  <-  pare  <-  FALSE
  phyp    <-  expand.grid('a'=rep(c(-1,1), each=3)*quantile(dat$vo2, probs=rep(c(0.025, 0.5, 0.975)), na.rm=TRUE), 'b'=rep(c(-1, 1), each=3)*rep(c(0.1,0.5,0.9), 2), 'c'=quantile(dat$vo2, probs=c(0.025, 0.5, 0.975), na.rm=TRUE))
  pmic    <-  ppar  <-  unique(phyp[,-3])

  n  <-  1
  while(!micmen) {
    cat('michaelisMenten\n')
    micmod  <-  try(nls(vo2 ~ michaelisMenten(a, b, po2), dat, start=pmic[n,]), silent=TRUE)
    if(class(micmod) == 'nls') {
      micmen  <-  TRUE
    } else {
      n  <-  n + 1
      if(n > nrow(pmic)) {
        micmod  <-  NA
        micmen  <-  TRUE
      }
    }
  }

  powmod  <-  lm(log(vo2) ~ log(po2), dat)

  n  <-  1
  while(!hyper) {
    cat('hyperbola\n')
    hypmod  <-  try(nls(vo2 ~ hyperbola(a, b, c, po2), dat, start=phyp[n,]), silent=TRUE)
    if(class(hypmod) == 'nls') {
      hyper  <-  TRUE
    } else {
      n  <-  n + 1
      if(n > nrow(phyp)) {
        hypmod  <-  NA
        hyper   <-  TRUE
      }
    }
  }
  n  <-  1
  while(!pare) {
    cat('pareto\n')
    parmod  <-  try(nls(vo2 ~ pareto(a, b, po2), dat, start=ppar[n,]), silent=TRUE)
    if(class(parmod) == 'nls') {
      pare  <-  TRUE
    } else {
      n  <-  n + 1
      if(n > nrow(ppar)) {
        parmod  <-  NA
        pare    <-  TRUE
      }
    }
  }

  list('michaelisMenten' = micmod,
       'powerFct' = powmod,
       'hyperbola' = hypmod,
       'pareto' = parmod)
}

compareAllModels  <-  function(allModels) {
  sapply(allModels, AIC)
}

calcPco2  <-  function(model, name) {
  coefs  <-  summary(model)$coefficients[1:2, 'Estimate']
  as.numeric(get(paste0('pco2', name))(coefs[1], coefs[2]))
}

pco2  <-  function(index, ...) {
  allModels  <-  runAllModels(...)
  allModels  <-  allModels[!is.na(allModels)]
  data.frame(column=index, fct=names(allModels), AICs=compareAllModels(allModels), pco2=mapply(calcPco2, allModels, names(allModels)), stringsAsFactors=FALSE, row.names=NULL)
}

pco2michaelisMenten  <-  function(a, b, m=0.065) {
  sqrt(a*b/m) - b
}

pco2powerFct  <-  function(a, b, m=0.065) {
  (m / (exp(a)*b)) ^ (1 / (b - 1))
}

pco2hyperbola  <-  function(a, b, m=0.065) {
  sqrt(a*b/m) - b
}

pco2pareto  <-  function(a, b, m=0.065) {
  (m / (-b * a^b)) ^ (1 / (b - 1))
}

