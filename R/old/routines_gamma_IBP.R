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
  
  return(list("means" = means,"ubs" = ub,"lbs" = lb))
  
}

#########################################################################

#' EB based on EFPF-max - IBP with Gamma mixture
#'
#' This function returns the value of (alpha, theta, a, b) maximizing the 
#' EFPF for the given sample - IBP with Gamma mixture
#'
#' @param n [integer] dimension of the observed sample
#' @param counts [numeric] vector of cardinalities for the observed features
#' @param pars_0 [numeric] Initialization for (alpha, theta, a, b) to 
#' be optimized
#' 
#' @export
EB_EFPF_gamma_IBP <- function(n, counts, pars_0){
  
  # initialize s = theta + alpha
  pars_0[2] <- pars_0[2] + pars_0[1]
  
  # set constraints stricter so that function is limited
  res <- optim(par = pars_0, fn = neg_log_EFPF_gamma_IBP_rep, n = n, counts = counts,
               method = "L-BFGS-B", lower = c(0,0.1,0.1, 0.1), 
               upper = c(0.99, Inf, Inf, Inf))
  
  sol <- res$par
  # convert to the theta parameter
  sol[2] <- sol[2] - sol[1]
  
  return (sol)
  
}
