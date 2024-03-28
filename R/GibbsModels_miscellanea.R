#
# Functions for predicting Kn + K^(n)_m, for given n, Kn and a range of M, a posteriori
# (i.e., using the posterior chains of the parameters).
# Available for both mixtures of Beta-Bernoulli and Gamma IBP.
#

#' Model-based prediction for BB with Poisson(lambda) mixture
#'
#' @param object An object of class \code{GibbsFA, PoissonBB}
#' @param n Size of the hypothetical observed sample
#' @param Kn Number of distinct features hypothetically observed in the sample
#' @param M Maximum size of future sample to extrapolate
#' @param seed Set seed for randomness
#' 
#' @details Draw samples from the posterior distribution Kn + K^(n)_m, for different m,
#' for BB with Poisson(lambda) mixture
prediction_PoissonBB <- function(object, n = 0, Kn = 0, M = 10, seed = 1234) {
  
  set.seed(seed)
  
  # Extract the chains and the parameters from object
  alpha_chain <- object$alpha_chain
  theta_chain <- object$theta_chain
  lambda <- object$prior$lambda
  
  # Number of iterations of the MCMC ( = number of samples from the posterior)
  S <- length(alpha_chain)
  
  M_vec <- 1:M
  
  # Draw samples from the posteriors and store in prediction_chain
  kmn_chain <- matrix(NA, nrow = M, ncol = S )
  
  for (q in 1:S){
    alpha <- alpha_chain[q]
    theta <- theta_chain[q]
    
    par_1 <- lambda*exp(lgamma(theta+alpha+n) - lgamma(theta+alpha) - 
                          lgamma(theta+n) + lgamma(theta))
    par_2 <- lambda*exp(lgamma(theta+alpha+n+M_vec) - lgamma(theta+alpha+n) - 
                          lgamma(theta+n+M_vec) + lgamma(theta+n) +
                          lgamma(theta+alpha+n) - lgamma(theta+alpha) - 
                          lgamma(theta+n) + lgamma(theta) )
    poiss_par <- par_1 - par_2
    
    kmn_chain[,q] <- rpois(M, poiss_par)
  }
  
  prediction_chain <- kmn_chain + Kn
  
  # Return a list 
  prediction_list <- as.list(data.frame(t(prediction_chain)))
  names(prediction_list) <- paste0("M = ", M_vec)
  
  return (prediction_list)
}


#' Model-based prediction for BB with NB(n0,mu0) mixture
#'
#' @param object An object of class \code{GibbsFA, NegBinBB}
#' @param n Size of the hypothetical observed sample
#' @param Kn Number of distinct features hypothetically observed in the sample
#' @param M Maximum size of future sample to extrapolate
#' @param seed Set seed for randomness
#' 
#' @details Draw samples from the posterior distribution Kn + K^(n)_m, for different m,
#' for BB with NB(n0,mu0) mixture
prediction_NegBinBB <- function(object, n = 0, Kn = 0, M = 10, seed = 1234) {
  
  set.seed(seed)
  # Extract the chains and the parameters from object
  alpha_chain <- object$alpha_chain
  theta_chain <- object$theta_chain
  n0 <- object$prior$n0
  mu0 <- object$prior$mu0
  
  # Number of iterations of the MCMC ( = number of samples from the posterior)
  S <- length(alpha_chain)
  
  M_vec <- 1:M
  
  # Draw samples from the posteriors and store in prediction_chain
  kmn_chain <- matrix(NA, nrow = M, ncol = S )
  
  for (q in 1:S){
    
    alpha <- alpha_chain[q]
    theta <- theta_chain[q]
    
    par_0 = (1-mu0)*exp(lgamma(theta) - lgamma(theta+alpha))
    par_1 = exp(lgamma( theta + alpha +n)- lgamma(theta +n))
    par_2 = exp(lgamma(theta+alpha+n+M_vec) - lgamma(theta+n+M_vec))
    
    p_bar <- 1- par_0 *(par_1 - par_2)/(1-par_0*par_2)
    
    kmn_chain[,q] <- rnbinom(M, n0 + Kn, p_bar)
  }
  
  prediction_chain <- kmn_chain + Kn
  
  # Return a list 
  prediction_list <- as.list(data.frame(t(prediction_chain)))
  names(prediction_list) <- paste0("M = ", M_vec)
  
  return (prediction_list)
}



#' Model-based prediction for IBP with Gamma(a,b) mixture,
#' with prior on a and b
#'
#' @param object An object of class \code{GibbsFA, GammaIBP}
#' @param n Size of the hypothetical observed sample
#' @param Kn Number of distinct features hypothetically observed in the sample
#' @param M Maximum size of future sample to extrapolate
#' @param seed Set seed for randomness
#' 
#' @details Draw samples from the posterior distribution Kn + K^(n)_m, for different m,
#' for IBP with Gamma(a,b) mixture, with prior on a and b
prediction_GammaIBP <- function(object, n = 0, Kn = 0, M = 10, seed = 1234) {
  
  set.seed(seed)
  
  # Extract the chains and the parameters from object
  alpha_chain <- object$alpha_chain
  theta_chain <- object$theta_chain
  a_chain <- object$a_chain
  b_chain <- object$b_chain
  
  # Number of iterations of the MCMC ( = number of samples from the posterior)
  S <- length(a_chain)
  
  M_vec <- 1:M
  
  # Draw samples from the posteriors and store in prediction_chain
  kmn_chain <- matrix(NA, nrow = M, ncol = S )
  
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
    
    sum_M <- stable_sum_M_all_gamma_IBP(alpha, theta, M, n)
    
    kmn_chain[,w] <- rnbinom(M, a + Kn, (gamma_a_t_n + b)/(gamma_a_t_n + b + sum_M))
    
    print(paste0("iteration: ", w))
    
  }
  
  prediction_chain <- kmn_chain + Kn
  
  # Return a list 
  prediction_list <- as.list(data.frame(t(prediction_chain)))
  names(prediction_list) <- paste0("M = ", M_vec)
  
  return (prediction_list)
}
