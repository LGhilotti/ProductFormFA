#
# Functions for predicting Kn + K^(n)_m, for given n, Kn and a range of M, a posteriori
# (i.e., using the posterior chains of the parameters).
# Available for both mixtures of Beta-Bernoulli and Gamma IBP.
# For the EB versions, the core parameters of the posterior distibutions are returned.
#

#' Model-based extrapolation
#' @param object An object of class \code{GibbsFA}
#' @param ... Additional parameters
#' 
#' @export
prediction <- function(object, ...) {
  
  UseMethod("prediction", object)
}


#' Model-based prediction for BB with Poisson(lambda) mixture
#'
#' @param object An object of class \code{GibbsFA, PoissonBB}
#' @param n Size of the hypothetical observed sample
#' @param Kn Number of distinct features hypothetically observed in the sample
#' @param M Maximum size of future sample to extrapolate
#' @param only_last logical; TRUE if only extrapolation for the additional sample of size M
#' @param seed Set seed for randomness
#' 
#' @details Draw samples from the posterior distribution Kn + K^(n)_m, for different m,
#' for BB with Poisson(lambda) mixture
prediction.PoissonBB <- function(object, n = 0, Kn = 0, M = 10, only_last = FALSE, seed = 1234) {
  
  set.seed(seed)
  
  # Extract the chains and the parameters from object
  alpha_chain <- object$alpha_chain
  theta_chain <- object$theta_chain
  lambda <- object$prior$lambda
  
  # Number of iterations of the MCMC ( = number of samples from the posterior)
  S <- length(alpha_chain)
  
  if (only_last == TRUE){
    M_vec <- M
  } else {
    M_vec <- 1:M
  }
  
  # Draw samples from the posteriors and store in prediction_chain
  kmn_chain <- matrix(NA, nrow = length(M_vec), ncol = S )
  
  for (q in 1:S){
    alpha <- alpha_chain[q]
    theta <- theta_chain[q]
    
    poiss_par <- lambda*feature_fraction(M_vec, alpha, theta + n)*
      (1 - feature_fraction(n, alpha, theta))
    
    kmn_chain[,q] <- rpois(length(M_vec), poiss_par)
  }
  
  prediction_chain <- kmn_chain + Kn
  
  # Return a list 
  prediction_list <- as.list(data.frame(t(prediction_chain)))
  names(prediction_list) <- paste0("M = ", M_vec)
  
  return (prediction_list)
}


#' Model-based prediction for BB with Poisson(lambda) mixture - EB version
#'
#' @param object An object of class \code{GibbsFA, PoissonBB}
#' @param n Size of the hypothetical observed sample
#' @param Kn Number of distinct features hypothetically observed in the sample
#' @param M Maximum size of future sample to extrapolate
#' @param only_last logical; TRUE if only extrapolation for the additional sample of size M
#' 
#' @details Return lambdas sequence of the posterior distribution K^(n)_m, for different m,
#' for BB with Poisson(lambda) mixture
prediction.PoissonBB_eb <- function(object, n = 0, Kn = 0, M = 10, only_last = FALSE) {
  
  # Extract the parameters from object
  alpha <- object$alpha
  theta <- object$theta
  lambda <- object$lambda

  if (only_last == TRUE){
    M_vec <- M
  } else {
    M_vec <- 1:M
  }
  
  poiss_par <- lambda*feature_fraction(M_vec, alpha, theta + n)*
    (1 - feature_fraction(n, alpha, theta))
    
  return (list("lambda_post" = poiss_par))
}



#' Model-based prediction for BB with NB(n0,mu0) mixture
#'
#' @param object An object of class \code{GibbsFA, NegBinBB}
#' @param n Size of the hypothetical observed sample
#' @param Kn Number of distinct features hypothetically observed in the sample
#' @param M Maximum size of future sample to extrapolate
#' @param only_last logical; TRUE if only extrapolation for the additional sample of size M
#' @param seed Set seed for randomness
#' 
#' @details Draw samples from the posterior distribution Kn + K^(n)_m, for different m,
#' for BB with NB(n0,mu0) mixture
prediction.NegBinBB <- function(object, n = 0, Kn = 0, M = 10, only_last = FALSE, seed = 1234) {
  
  set.seed(seed)
  # Extract the chains and the parameters from object
  alpha_chain <- object$alpha_chain
  theta_chain <- object$theta_chain
  n0 <- object$prior$n0
  mu0 <- object$prior$mu0
  
  # Number of iterations of the MCMC ( = number of samples from the posterior)
  S <- length(alpha_chain)
  
  if (only_last == TRUE){
    M_vec <- M
  } else {
    M_vec <- 1:M
  }
  
  # Draw samples from the posteriors and store in prediction_chain
  kmn_chain <- matrix(NA, nrow = length(M_vec), ncol = S )
  
  for (q in 1:S){
    
    alpha <- alpha_chain[q]
    theta <- theta_chain[q]
    
    # in (nstar,p) parametrization
    n0_upd <- n0 + Kn
    mu0_upd <- (n0 + Kn)/(n0/mu0 + feature_fraction(n, alpha, theta) ) *
      feature_fraction(M_vec, alpha, theta + n) * (1 - feature_fraction(n, alpha, theta))
    
    p_upd <- 1/(mu0_upd/n0_upd + 1)
    
    # par_0 = (1-p)*exp(lgamma(theta) - lgamma(theta+alpha))
    # par_1 = exp(lgamma( theta + alpha +n)- lgamma(theta +n))
    # par_2 = exp(lgamma(theta+alpha+n+M_vec) - lgamma(theta+n+M_vec))
    # 
    # p_bar <- 1- par_0 *(par_1 - par_2)/(1-par_0*par_2)
     
    kmn_chain[,q] <- rnbinom(length(M_vec), n0_upd, p_upd)
  }
  
  prediction_chain <- kmn_chain + Kn
  
  # Return a list 
  prediction_list <- as.list(data.frame(t(prediction_chain)))
  names(prediction_list) <- paste0("M = ", M_vec)
  
  return (prediction_list)
}


#' Model-based prediction for BB with NB(n0,mu0) mixture - EB version
#'
#' @param object An object of class \code{GibbsFA, NegBinBB}
#' @param n Size of the hypothetical observed sample
#' @param Kn Number of distinct features hypothetically observed in the sample
#' @param M Maximum size of future sample to extrapolate
#' @param only_last logical; TRUE if only extrapolation for the additional sample of size M
#' 
#' @details Return the (n0,mu0)s sequences of the NB of K^(n)_m, for different m,
#' for BB with NB(n0,mu0) mixture
prediction.NegBinBB_eb <- function(object, n = 0, Kn = 0, M = 10, only_last = FALSE) {
  
  # Extract the parameters from object
  alpha <- object$alpha
  theta <- object$theta
  n0 <- object$n0
  mu0 <- object$mu0
  
  p_n <- feature_fraction(n, alpha, theta)
  
  if (only_last == TRUE){
    M_vec <- M
  } else {
    M_vec <- 1:M
  }
  
  n0_upd <- n0 + Kn
  mu0_upd <- (n0 + Kn)/(n0/mu0 + feature_fraction(n, alpha, theta) ) *
    feature_fraction(M_vec, alpha, theta + n) * (1 - feature_fraction(n, alpha, theta))
  
  
  return (list("n0_post" = rep(n0_upd, length(M_vec)),
               "mu0_post" = mu0_upd) )
}




#' Model-based prediction for IBP with Gamma(a,b) mixture,
#' with prior on a and b
#'
#' @param object An object of class \code{GibbsFA, GammaIBP}
#' @param n Size of the hypothetical observed sample
#' @param Kn Number of distinct features hypothetically observed in the sample
#' @param M Maximum size of future sample to extrapolate
#' @param only_last logical; TRUE if only extrapolation for the additional sample of size M
#' @param seed Set seed for randomness
#' 
#' @details Draw samples from the posterior distribution Kn + K^(n)_m, for different m,
#' for IBP with Gamma(a,b) mixture, with prior on a and b
prediction.GammaIBP <- function(object, n = 0, Kn = 0, M = 10, only_last = FALSE, seed = 1234) {
  
  set.seed(seed)
  
  # Extract the chains and the parameters from object
  alpha_chain <- object$alpha_chain
  theta_chain <- object$theta_chain
  a_chain <- object$a_chain
  b_chain <- object$b_chain
  
  # Number of iterations of the MCMC ( = number of samples from the posterior)
  S <- length(a_chain)
  
  if (only_last == TRUE){
    M_vec <- M
  } else {
    M_vec <- 1:M
  }
  
  # Draw samples from the posteriors and store in prediction_chain
  kmn_chain <- matrix(NA, nrow = length(M_vec), ncol = S )
  
  for (w in 1:S){
    a <- a_chain[w]
    b <- b_chain[w]
    alpha <- alpha_chain[w]
    theta <- theta_chain[w]
    
    if (n == 0){
      gamma_a_t_n <- 0
    } else {
      gamma_a_t_n <- sum(exp(lgamma(alpha + theta + (1:n) - 1) - lgamma(alpha + theta) -
                               lgamma(theta + (1:n)) + lgamma(theta +1) ) )
    }
    
    sum_M <- stable_sum_M_all_gamma_IBP(alpha, theta, M, only_last, n)
    
    kmn_chain[,w] <- rnbinom(length(M_vec), a + Kn, (gamma_a_t_n + b)/(gamma_a_t_n + b + sum_M))
    
    print(paste0("iteration: ", w))
    
  }
  
  prediction_chain <- kmn_chain + Kn
  
  # Return a list 
  prediction_list <- as.list(data.frame(t(prediction_chain)))
  names(prediction_list) <- paste0("M = ", M_vec)
  
  return (prediction_list)
}


#' Model-based prediction for IBP with Gamma(a,b) mixture - EB version
#'
#' @param object An object of class \code{GibbsFA, GammaIBP}
#' @param n Size of the hypothetical observed sample
#' @param Kn Number of distinct features hypothetically observed in the sample
#' @param M Maximum size of future sample to extrapolate
#' @param only_last logical; TRUE if only extrapolation for the additional sample of size M
#' 
#' @details Return the (n0,mu0) sequences of the NB posterior distribution K^(n)_m, 
#' for different m, for IBP with Gamma(a,b) mixture
prediction.GammaIBP_eb <- function(object, n = 0, Kn = 0, M = 10, only_last = FALSE) {
  
  # Extract the parameters from object
  alpha <- object$alpha
  theta <- object$theta
  a <- object$a
  b <- object$b
  
  if (only_last == TRUE){
    M_vec <- M
  } else {
    M_vec <- 1:M
  }
  
  if (n == 0){
    gamma_a_t_n <- 0
  } else {
    gamma_a_t_n <- sum(exp(lgamma(alpha + theta + (1:n) - 1) - lgamma(alpha + theta) -
                             lgamma(theta + (1:n)) + lgamma(theta +1) ) )
  }
  
  sum_M <- stable_sum_M_all_gamma_IBP(alpha, theta, M, only_last, n)
    
  n0_upd <- a + Kn
  p_upd <- (gamma_a_t_n + b)/(gamma_a_t_n + b + sum_M)
  mu0_upd <- n0_upd*(1 - p_upd)/p_upd
    
    
  return (list("n0_post" = rep(n0_upd, length(M_vec)),
               "mu0_post" = mu0_upd))
}