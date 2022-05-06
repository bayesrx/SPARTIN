#' Samples from true posterior of modified Hierarchical Strauss,
#' as well as simulated null distribution. Used in Computation of CTIP.
#' @param tile ppp object with two qualitative marks
#' @param r radius of interaction for model fitting.
#' @param quad.spacing space between points in quadrature used in estimation.
#' @param correction correction used in fitting; see spatstat documentation of
#' "ppm.ppp" for more details.
#' @param n.null number of null simulations; larger numbers of null simulations
#' improve precision, but are computationally more demanding
#' @param n.burn number of burn-in samples for MCMC model fitting for actual
#' tile as well as null tiles
#' @param n.sample total number of samples drawn for each simulation. The total
#' number of samples for the actual tile will be (n.sample - n.burn)/n.thin
#' @param n.thin thinning for MCMC fitting
#' @param log.beta.1.mean prior mean for log(beta_1), i.e. log first order
#' intensity of type 1 points
#' @param log.beta.1.prec prior precision for log(beta_1)
#' @param log.beta.2.mean prior mean for log(beta_2), i.e. log first order
#' intensity of type 2 points
#' @param log.beta.2.prec prior precision for log(beta_2)
#' @param log.gamma.mean prior mean for log(gamma), i.e. interaction parameter
#' @param log.gamma.prec prior precision for log(gamma)
#' @return list with three attributes: "freq" gives a summary of the frequentist
#' fitting (for diagnostics), "model" gives the R2Jags object corresponding
#' to the actual fitted model, and "null.samples" gives samples from simulated
#' null distribution
#' @noRd
FullFit = function(tile, r, quad.spacing,
                   n.null = 5,
                   n.burn = 1000, n.sample = 11000,
                   null.n.burn = 1000, null.n.sample = 11000,
                   n.thin = 5,
                   log.beta.1.mean = 0, log.beta.1.prec = 0.0000001,
                   log.beta.2.mean = 0,  log.beta.2.prec = 0.0000001,
                   log.gamma.mean = 0, log.gamma.prec = 0.0000001,
                   correction = "Ripley"){

  freq.model = FitHSFreq(tile, r = r, quad.spacing = quad.spacing,
                         correction = correction)
  #
  # glmdata = freq.model$internal$glmdata

  window.area = area.owin(tile$window)
  tumor.intensity = sum(tile$marks == 1)/window.area
  til.intensity = sum(tile$marks == 2)/window.area

  tums = ppp(x = tile$x[tile$marks == 1],
             y = tile$y[tile$marks == 1],
             window = tile$window)



  bayesian.model = FitHSBayes(tile,
                              r = r, quad.spacing = quad.spacing,
                              n.burn = n.burn, n.sample = n.sample,
                              n.thin = n.thin,
                              log.beta.1.mean = log.beta.1.mean,
                              log.beta.1.prec = log.beta.1.prec,
                              log.beta.2.mean = log.beta.2.mean,
                              log.beta.2.prec = log.beta.2.prec,
                              log.gamma.mean = log.gamma.mean,
                              log.gamma.prec = log.gamma.prec)

  null.samples = unlist(map(1:n.null, function(x){
    tils = rpoispp(til.intensity, win = tile$window)

    x = c(tums$x, tils$x)
    y = c(tums$y, tils$y)
    marks = factor(c(rep(1, tums$n), rep(2, tils$n)))

    null.sim = ppp(x = x, y = y, marks = marks, window = tile$window)
    # stop()
    # null.freq = FitHSFreq(null.sim, r = r, quad.spacing = quad.spacing,
    #                       correction = correction)
    null.model = FitHSBayes(null.sim, # null.freq$internal$glmdata,
                            r = r, quad.spacing = quad.spacing,
                            n.burn = null.n.burn, n.sample = null.n.sample,
                            n.thin = n.thin,
                            log.beta.1.mean = log.beta.1.mean,
                            log.beta.1.prec = log.beta.1.prec,
                            log.beta.2.mean = log.beta.2.mean,
                            log.beta.2.prec = log.beta.2.prec,
                            log.gamma.mean = log.gamma.mean,
                            log.gamma.prec = log.gamma.prec)

    return(null.model$BUGSoutput$sims.list$log.gamma)
  }))

  ret.obj = {}

  ret.obj$freq = summary(freq.model)$coefs.SE.CI
  ret.obj$model = bayesian.model
  ret.obj$null.samples = null.samples

  return(ret.obj)
}
