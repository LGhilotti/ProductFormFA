#' Credible intervals of Kmn for the IBP with Gamma mixture
#'
#' This function computes the means and the credible intervals of Kmn for the IBP with Gamma mixture
#'
#' @param alpha [numeric] value of alpha in product-form feature allocation
#' @param theta [numeric] value of theta in product-form feature allocation
#' @param m [numeric] dimension of the new sample to be observed
#' @param n [numeric] dimension of the already observed sample
#' @param Kn [numeric] number of features in the already observed sample
#' @param a [numeric] Gamma hyperparameter (shape)
#' @param b [numeric] Gamma hyperparameter (rate)
#' @param lev [numeric] level of the credible intervals
#'
#' @export
#'
CI_Kmn_gamma_IBP <- function(alpha, theta, m, n, Kn, a, b, lev) {
  
  pbars <- p_kmn_all_gamma_IBP(alpha, theta, m, n, b)
  means <- (a + Kn) * (1 - pbars) / pbars
  ub <- qnbinom(lev + (1 - lev) / 2, a + Kn, pbars)
  lb <- qnbinom((1 - lev) / 2, a + Kn, pbars)
  
  return(list(means,ub,lb))
  
}
