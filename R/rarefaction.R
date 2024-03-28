#
# Functions for computing the rarefaction curve, i.e. the statistic K_n, a posteriori
# (i.e., using the posterior chains of the parameters).
# Available for both mixtures of Beta-Bernoulli and Gamma IBP.
#

#' Model-based rarefaction
#' @param object An object of class \code{GibbsFA}
#' @param ... Additional parameters
#' 
#' @export
rarefaction <- function(object, ...) {
  
  UseMethod("rarefaction", object)
}


#' Sample-based rarefaction curve for a vector of species abundances 
#' (RIGHT NOW IT IS JUST THE OBSERVED SAMPLE, WITH NO AVERAGE)
#' 
#' @param object  The \code{n x K}-dimensional binary matrix
#' 
#' @export
#' 
rarefaction.array <- function(object) {
  
  feature_list <- convert_features_list(object)
  n <- nrow(object)
  
  Kn_observed <- sapply(1:n, function(i) length(unique(unlist(feature_list[1:i]))) )
  
  return(Kn_observed)
  
}


#' Model-based rarefaction for BB with Poisson(lambda) mixture
#'
#' @param object An object of class \code{GibbsFA, PoissonBB}
#' @param seed Set seed for randomness
#' 
#' @export
#' 
#' @details Draw samples from the posterior distribution of K_n,
#' for BB with Poisson(lambda) mixture
rarefaction.PoissonBB <- function(object, seed = 1234) {
  
  # Extract the size of the observed sample "n"
  feature_matrix <- object$feature_matrix
  n <- nrow(feature_matrix)
  
  # Perform the prediction
  rarefaction_list <- prediction_PoissonBB(object, n = 0, Kn = 0, M = n, seed = seed)
  
  names(rarefaction_list) <- paste0("N = ", 1:n)
  
  return (rarefaction_list)
}


#' Model-based rarefaction for BB with NB(n0,mu0) mixture
#'
#' @param object An object of class \code{GibbsFA, NegBinBB}
#' @param seed Set seed for randomness
#' 
#' @export
#' 
#' @details Draw samples from the posterior distribution of K_n,
#' for BB with NB(n0,mu0) mixture
rarefaction.NegBinBB <- function(object, seed = 1234) {
  
  # Extract the size of the observed sample "n"
  feature_matrix <- object$feature_matrix
  n <- nrow(feature_matrix)

  # Perform the prediction
  rarefaction_list <- prediction_NegBinBB(object, n = 0, Kn = 0, M = n, seed = seed)
  
  names(rarefaction_list) <- paste0("N = ", 1:n)
  
  return (rarefaction_list)
}




#' Model-based rarefaction for IBP with Gamma(a,b) mixture,
#' with prior on a and b
#'
#' @param object An object of class \code{GibbsFA, GammaIBP}
#' @param seed Set seed for randomness
#' 
#' @export
#' 
#' @details Draw samples from the posterior distribution of K_n,
#' for IBP with Gamma(a,b) mixture, and prior on a and b
rarefaction.GammaIBP <- function(object, seed = 1234) {
  
  # Extract the size of the observed sample "n"
  feature_matrix <- object$feature_matrix
  n <- nrow(feature_matrix)
  
  # Perform the prediction
  rarefaction_list <- prediction_GammaIBP(object, n = 0, Kn = 0, M = n, seed = seed)
  
  names(rarefaction_list) <- paste0("N = ", 1:n)
  
  return (rarefaction_list)
}
