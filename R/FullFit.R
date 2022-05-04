#' Samples from posterior

FullFit = function(tile, n.null = 5, n.burn = 1000, n.sample = 11000,
                   n.thin = 5,
                   log.beta.1.mean = 0, log.beta.1.prec = 0.0000001,
                   log.beta.2.mean = 0,  log.beta.2.prec = 0.0000001,
                   log.gamma.mean = 0, log.gamma.prec = 0.0000001){
  freq.model = FitHSPPPFull(tile)
  glmdata = freq.model$internal$glmdata

  window.area = area.owin(tile$window)
  tumor.intensity = sum(tile$marks == 1)/window.area
  til.intensity = sum(tile$marks == 2)/window.area

  tums = ppp(x = tile$x[tile$marks == 1],
             y = tile$y[tile$marks == 1],
             window = tile$window)



  bayesian.model = PBHSFast(glmdata)

  null.samples = unlist(map(1:n.null, function(x){
    tils = rpoispp(til.intensity, win = tile$window)

    x = c(tums$x, tils$x)
    y = c(tums$y, tils$y)
    marks = factor(c(rep(1, tums$n), rep(2, tils$n)))

    null.sim = ppp(x = x, y = y, marks = marks, window = tile$window)
    # stop()
    null.freq = FitHSPPPFull(null.sim, r = 30)
    null.model = PBHSFast(null.freq$internal$glmdata)

    return(null.model$BUGSoutput$sims.list$log.gamma,
           n.burn = n.burn, n.sample = n.sample,
           n.thin = n.thin,
           log.beta.1.mean = log.beta.1.mean,
           log.beta.1.prec = log.beta.1.prec,
           log.beta.2.mean = log.beta.2.mean,
           log.beta.2.prec = log.beta.2.prec,
           log.gamma.mean = log.gamma.mean,
           log.gamma.prec = log.gamma.prec)
  }))

  ret.obj = {}

  ret.obj$freq = summary(freq.model)$coefs.SE.CI
  ret.obj$model = bayesian.model
  ret.obj$null.samples = null.samples

  return(ret.obj)
}
