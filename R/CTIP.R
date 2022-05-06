#' Computes CTIP for a given tile
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
#' number of samples for the actual tile will be (n.sample - n.burn)/n.thin.
#' This also determines how many samples will be used to perform stochastic
#' integration used to estimate CTIP
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

CTIP = function(tile, r, quad.spacing,
                correction = "Ripley",
                n.null = 5,
                n.burn = 1000, n.sample = 11000,
                null.n.burn = 1000, null.n.sample = 11000,
                n.thin = 5,
                log.beta.1.mean = 0, log.beta.1.prec = 0.0000001,
                log.beta.2.mean = 0,  log.beta.2.prec = 0.0000001,
                log.gamma.mean = 0, log.gamma.prec = 0.0000001){

  ff = FullFit(tile,
               r = r, quad.spacing = quad.spacing,
               correction = correction,
               n.null = n.null,
               n.burn = n.burn, n.sample = n.sample,
               null.n.burn = null.n.burn, null.n.sample = null.n.sample,
               n.thin = n.thin,
               log.beta.1.mean = log.beta.1.mean,
               log.beta.1.prec = log.beta.1.prec,
               log.beta.2.mean = log.beta.2.mean,
               log.beta.2.prec = log.beta.2.prec,
               log.gamma.mean = log.gamma.mean,
               log.gamma.prec = log.gamma.prec)

  actual_samples = ff$model$BUGSoutput$sims.list$log.gamma
  if(length(ff$null.samples) < length(actual_samples)){
    null_samples = sample(ff$null.samples,
                          length(actual_samples),
                          replace = TRUE)
  }else{
    null_samples = sample(ff$null.samples,
                          length(actual_samples))
  }

  return(mean(actual_samples > null_samples))
}
