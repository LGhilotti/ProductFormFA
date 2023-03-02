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
