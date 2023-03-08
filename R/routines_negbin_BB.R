#' Credible intervals of Kmn for the BB with Negative-Binomial mixture
#'
#' This function computes the means and the credible intervals of Kmn for the BB with Negative-Binomial mixture
#'
#' @param alpha [numeric] value of alpha in product-form feature allocation
#' @param theta [numeric] value of theta in product-form feature allocation
#' @param m [numeric] dimension of the new sample to be observed
#' @param n [numeric] dimension of the already observed sample
#' @param Kn [numeric] number of features in the already observed sample
#' @param nstar [numeric] Negative-Binomial hyperparameter (number of successes)
#' @param p [numeric] Negative-Binomial hyperparameter (success probability)
#' @param lev [numeric] level of the credible intervals
#'
#' @export
#'
CI_Kmn_negbin_BB <- function(alpha, theta, m, n, Kn, nstar, p, lev) {
  
  pbars <- p_kmn_all_negbin_BB(alpha, theta, m, n, p)

  means <- (nstar + Kn) * (1 - pbars) / pbars

  ub <- qnbinom(lev + (1 - lev) / 2, nstar + Kn, pbars)

  lb <- qnbinom((1 - lev) / 2, nstar + Kn, pbars)

  return(list("means" = means,"ubs" = ub,"lbs" = lb))
  
}

#########################################################################

#' EB based on EFPF-max - BB with Negative-Binomial mixture
#'
#' This function returns the value of (alpha, theta, n*, p) maximizing the 
#' EFPF for the given sample - BB with Negative-Binomial mixture
#'
#' @param n [integer] dimension of the observed sample
#' @param counts [numeric] vector of cardinalities for the observed features
#' @param pars_0 [numeric] Initialization for (alpha, theta, n*, p) to 
#' be optimized
#' 
#' @export
EB_EFPF_negbin_BB <- function(n, counts, pars_0){
  
  # initialize s = theta + alpha
  pars_0[2] <- pars_0[2] + pars_0[1]
  
  # set constraints stricter so that function is limited
  res <- optim(par = pars_0, fn = neg_log_EFPF_negbin_BB_rep, n = n, counts = counts,
               method = "L-BFGS-B", lower = c(-Inf,0.1,0.1, 0.01), 
               upper = c(-0.1, Inf, Inf, 0.99))
  
  sol <- res$par
  # convert to the theta parameter
  sol[2] <- sol[2] - sol[1]
  
  return (sol)
  
}
