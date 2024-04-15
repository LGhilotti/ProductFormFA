#
# Functions for total richness estimation.
# Available only for mixtures of Beta-Bernoulli models.
#

#' Model-based richness posterior distribution
#' @param object An object of class \code{GibbsFA}
#' @param ... Additional parameters
#' 
#' @export
total_richness <- function(object, ...) {
  
  UseMethod("total_richness", object)
}



feature_fraction <- function(n, alpha, theta){
  1 - exp(lgamma(theta + alpha + n ) - lgamma(theta + alpha) - lgamma(theta + n) + lgamma(theta))
}

feature_fraction <- Vectorize(feature_fraction, "n")


#' Model-based richness posterior distribution for BB with Poisson(lambda) mixture
#'
#' @param object An object of class \code{GibbsFA, PoissonBB}
#' @param seed Set seed for randomness
#' 
#' @export
#' 
#' @details Draw samples from the posterior distribution of the total richness
#' for BB with Poisson(lambda) mixture
total_richness.PoissonBB <- function(object, seed = 1234) {
  
  set.seed(seed)
  
  # Extract the size of the observed sample "n" and the number of distinct features "Kn"
  feature_matrix <- object$feature_matrix
  n_train <- nrow(feature_matrix)
  Kn <- ncol(feature_matrix)
  
  # Extract the chains and the parameters from object
  alpha_chain <- object$alpha_chain
  theta_chain <- object$theta_chain
  lambda <- object$prior$lambda
  
  # Number of iterations of the MCMC ( = number of samples from the posterior)
  S <- length(alpha_chain)
  
  # Draw samples from the posterior and store in richness_chain
  richness_chain <- vector(length = S )
  
  for (q in 1:S){
    alpha <- alpha_chain[q]
    theta <- theta_chain[q]
    
    poiss_par <- lambda*(1 - feature_fraction(n_train, alpha, theta))
    
    richness_chain[q] <- Kn + rpois(1, poiss_par)
  }
  
  return (richness_chain)
}



#' Model-based richness posterior distribution for BB with Poisson(lambda) mixture - EB version
#'
#' @param object An object of class \code{GibbsFA, PoissonBB_eb}
#' 
#' @export
#' 
#' @details Return the lambda parameter of the Poisson posterior distribution of N'
#' for BB with Poisson(lambda) mixture
total_richness.PoissonBB_eb <- function(object) {
  
  alpha <- object$alpha
  theta <- object$theta
  lambda <- object$lambda
  n <- nrow(object$feature_matrix)
  
  return( list("lambda_post" = lambda * (1 - feature_fraction(n, alpha, theta)) ) )
}


#' Model-based richness posterior distribution for BB with NB(n0,mu0) mixture
#'
#' @param object An object of class \code{GibbsFA, NegBinBB}
#' @param seed Set seed for randomness
#' 
#' @export
#' 
#' @details Draw samples from the posterior distribution of the total richness
#' for BB with NB(n0,mu0) mixture
total_richness.NegBinBB <- function(object, seed = 1234) {
  
  set.seed(seed)
  
  # Extract the size of the observed sample "n" and the number of distinct features "Kn"
  feature_matrix <- object$feature_matrix
  n_train <- nrow(feature_matrix)
  Kn <- ncol(feature_matrix)
  
  # Extract the chains and the parameters from object
  alpha_chain <- object$alpha_chain
  theta_chain <- object$theta_chain
  n0 <- object$prior$n0
  mu0 <- object$prior$mu0
  
  # Number of iterations of the MCMC ( = number of samples from the posterior)
  S <- length(alpha_chain)
  
  # Draw samples from the posterior and store in richness_chain
  richness_chain <- vector(length = S )
  
  for (q in 1:S){
    
    alpha <- alpha_chain[q]
    theta <- theta_chain[q]
    
    # in (nstar,p) parametrization
    p <- 1/(mu0/n0 + 1)
    negbin_par <- 1 - (1 - p)*(1 - feature_fraction(n_train, alpha, theta))
    
    richness_chain[q] <- Kn + rnbinom(1, n0 + Kn, negbin_par)
  }
  
  return (richness_chain)
}

#' Model-based richness posterior distribution for BB with NB(n0,mu0) mixture - EB version
#'
#' @param object An object of class \code{GibbsFA, NegBinBB_eb}
#' 
#' @export
#' 
#' @details Return the (n0,mu0) parameter of the NB posterior distribution of N'
#' for BB with NB(n0,mu0) mixture
total_richness.NegBinBB_eb <- function(object) {
  
  alpha <- object$alpha
  theta <- object$theta
  n0 <- object$n0
  mu0 <- object$mu0
  n <- nrow(object$feature_matrix)
  Kn <- ncol(object$feature_matrix)
  p_n <- feature_fraction(n, alpha, theta)
  
  return( list("n0_post" = n0 + Kn, 
               "mu0_post" = (n0 + Kn)/(n0/mu0 + p_n) *(1 - p_n) ) )
}