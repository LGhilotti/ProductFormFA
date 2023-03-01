#' Credible intervals of Kmn for the BB with Poisson mixture
#'
#' This function computes the means and the credible intervals of Kmn for the BB with Poisson mixture
#'
#' @param alpha [numeric] value of alpha in product-form feature allocation
#' @param theta [numeric] value of theta in product-form feature allocation
#' @param m [numeric] dimension of the new sample to be observed
#' @param n [numeric] dimension of the already observed sample
#' @param lambda [numeric] Poisson hyperparameter 
#' @param lev [numeric] level of the credible intervals
#'
#' @export
#'
CI_Kmn_poiss_BB <- function(alpha, theta, m, n, lambda, lev) {
  
  means <- mean_kmn_all_poiss_BB(alpha, theta, m, n, lambda)
  ub <- qpois(lev + (1 - lev) / 2, means)
  lb <- qpois((1 - lev) / 2, means)

  return(list("means" = means,"ubs" = ub,"lbs" = lb))
}
