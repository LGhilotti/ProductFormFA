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

#########################################################################

#' EB based on EFPF-max - BB with Poisson mixture
#'
#' This function returns the value of (alpha, theta, lambda) maximizing the 
#' EFPF for the given sample - BB with Poisson mixture
#'
#' @param n [integer] dimension of the observed sample
#' @param counts [numeric] vector of cardinalities for the observed features
#' @param pars_0 [numeric] Initialization for (alpha, theta, lambda) to 
#' be optimized
#' 
#' @export
EB_EFPF_poiss_BB <- function(n, counts, pars_0){
  
  # initialize s = theta + alpha
  pars_0[2] <- pars_0[2] + pars_0[1]
  
  # set constraints stricter so that function is limited
  res <- optim(par = pars_0, fn = neg_log_EFPF_poiss_BB_rep, n = n, counts = counts,
               method = "L-BFGS-B", lower = c(-Inf,0.1,0.1), upper = c(-0.1, Inf, Inf))
  
  sol <- res$par
  # convert to the theta parameter
  sol[2] <- sol[2] - sol[1]
  
  return (sol)
  
}
