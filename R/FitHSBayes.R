#' Fits a bayesian version of a special case of the Hierarchical Strauss Model
#' with gamma_{11} = gamma_{22} = 1, i.e. assuming no
#' interaction between points of the same type.
#'
#' @param p spatstat ppp object with two qualitative mark levels, 1 and 2.
#' Note that model must have been fit using Berman-Turner Device for this to work.
#' @param r radius of interaction for model fitting.
#' @param quad.spacing space between points in quadrature used in estimation.
#' @param correction correction used in fitting; see spatstat documentation of
#' "ppm.ppp" for more details.
#' @param n.chains Number of chains to run
#' @param n.sample Total number of posterior samples
#' @param n.burn Number of burn-in samples
#' @param n.thin Thinning for chains, i.e. chains will keep sample every n.thin-th
#' sample from each model fitting.
#' @param log.beta.1.mean Mean of normal prior on log(beta_1), first order
#' intensity of tumor cells
#' @param log.beta.1.prec *Precision* of normal prior on log(beta_1), the
#' log-first order intensity of tumor cells. Note smaller values
#' indicate a less informative prior.
#' @param log.beta.2.mean Defined analogously to log.beta.1.mean, but for
#' log(beta_2), the log-first order intensity of lymphocytes.
#' @param log.beta.2.prec Defined analogously to log.beta.1.prec, but for
#' log(beta_2), the log-first order intensity of lymphocytes.
#' @param log.gamma.mean Mean of normal prior on log(gamma_{12}), the
#' log-interaction parameter between tumor cells and lymphocytes. Higher values
#' indicate more positive interaction, lower values indicate
#' more negative interaction.
#' @param log.gamma.prec *Precision* of normal prior on log(gamma_{12}), the
#' log-interaction parameter between tumor cells and lymphocytes. Note that
#' smaller values indicate a less informative prior.
#' @return R2jags object
#' @import R2jags
#' @export
FitHSBayes = function(p, r, quad.spacing, correction = "Ripley",
                    n.chains = 1,
                    n.sample = 11000, n.burn = 1000,
                    n.thin = 5,
                    log.beta.1.mean = 0, log.beta.1.prec = 0.0000001,
                    log.beta.2.mean = 0,  log.beta.2.prec = 0.0000001,
                    log.gamma.mean = 0, log.gamma.prec = 0.0000001){

  freq_mod = FitHSFreq(p, r = r, quad.spacing = quad.spacing,
                       correction = correction)

  glmdata = freq_mod$internal$glmdata

  bayesian.hierstrauss = function(){
    for(i in 1:N){
      const.vec[i] = exp(log.beta.1*beta.1.ind[i])*exp(log.beta.2*(1 - beta.1.ind[i]))*exp(log.gamma*nn[i])*w[i]
    }

    lik.const = (1/n)*sum(const.vec)

    for(i in 1:n){
      spy[i] = -1*(log.beta.1*beta.1.ind[i] + log.beta.2*(1 - beta.1.ind[i])
                   + log.gamma*nn[i]
                   - lik.const) + 1000

      zeros[i] ~ dpois(spy[i])
    }

    log.beta.1 ~ dnorm(log.beta.1.mean, log.beta.1.prec)
    log.beta.2 ~ dnorm(log.beta.2.mean, log.beta.2.prec)
    log.gamma ~ dnorm(log.gamma.mean, log.gamma.prec)
  }

  y = glmdata$.mpl.Y
  beta.1.ind = as.numeric(glmdata$marks == 1)
  w = glmdata$.mpl.W
  nn = glmdata$markX1xX2
  N = nrow(glmdata)
  n = max(which(glmdata$.mpl.Y != 0))
  zeros = rep(0, n)

  inits.JAGS = vector('list', n.chains)
  for(i in 1:length(inits.JAGS)){
    inits.JAGS[[i]] = list(log.beta.1=0,log.beta.2=0,log.gamma=0.0)
  }

  para.JAGS = c("log.beta.1", "log.beta.2", "log.gamma")

  dat.JAGS = list(
    beta.1.ind = beta.1.ind,
    w = w,
    nn = nn,
    N = N,
    n = n,
    zeros = zeros,
    log.gamma.mean = log.gamma.mean,
    log.gamma.prec = log.gamma.prec,
    log.beta.1.mean = log.beta.1.mean,
    log.beta.1.prec = log.beta.1.prec,
    log.beta.2.mean = log.beta.2.mean,
    log.beta.2.prec = log.beta.2.prec
  )

  capture.output({
    fit.JAGS = jags(data=dat.JAGS, inits=inits.JAGS,
                    parameters.to.save = para.JAGS, n.chains=n.chains,
                    n.iter=n.sample, n.thin = n.thin,
                    n.burnin=n.burn, model.file=bayesian.hierstrauss)
  })


  return(fit.JAGS)
}
