#
# Functions for computing the extrapolation curve, i.e. the statistic k_n + K^(n)_m, a posteriori
# (i.e, using the posterior chains of the parameters).
# Available for both mixtures of Beta-Bernoulli and Gamma IBP.
#

#' Model-based extrapolation
#' @param object An object of class \code{GibbsFA}
#' @param ... Additional parameters
#' 
#' @export
extrapolation <- function(object, ...) {
  
  UseMethod("extrapolation", object)
}




#' Model-based extrapolation for BB with Poisson(lambda) mixture
#'
#' @param object An object of class \code{GibbsFA, PoissonBB}
#' @param M Maximum size of future sample to extrapolate
#' @param only_last logical; TRUE if only extrapolation for the additional sample of size M
#' @param seed Set seed for randomness
#' 
#' @export
#' 
#' @details Draw samples from the posterior distribution of k_n + K^(n)_m, for different m,
#' for BB with Poisson(lambda) mixture
extrapolation.PoissonBB <- function(object, M = 10, only_last = FALSE, seed = 1234) {
  
  # Extract the size of the observed sample "n" and the number of distinct features "Kn"
  feature_matrix <- object$feature_matrix
  n <- nrow(feature_matrix)
  Kn <- ncol(feature_matrix)
  
  # Perform the prediction
  extrapolation_list <- prediction_PoissonBB(object, n = n, Kn = Kn, M = M, only_last = only_last, seed = seed)
  
  return (extrapolation_list)
}


#' Model-based extrapolation for BB with NB(n0,mu0) mixture
#'
#' @param object An object of class \code{GibbsFA, NegBinBB}
#' @param M Maximum size of future sample to extrapolate
#' @param only_last logical; TRUE if only extrapolation for the additional sample of size M
#' @param seed Set seed for randomness
#' 
#' @export
#' 
#' @details Draw samples from the posterior distribution of k_n + K^(n)_m, for different m,
#' for BB with NB(n0,mu0) mixture
extrapolation.NegBinBB <- function(object, M = 10, only_last = FALSE, seed = 1234) {
  
  # Extract the size of the observed sample "n" and the number of distinct features "Kn"
  feature_matrix <- object$feature_matrix
  n <- nrow(feature_matrix)
  Kn <- ncol(feature_matrix)
  
  # Perform the prediction
  extrapolation_list <- prediction_NegBinBB(object, n = n, Kn = Kn, M = M, only_last = only_last, seed = seed)
  
  return (extrapolation_list)
}





#' Model-based extrapolation for IBP with Gamma(a,b) mixture,
#' with prior on a and b
#'
#' @param object An object of class \code{GibbsFA, GammaIBP}
#' @param M Maximum size of future sample to extrapolate
#' @param only_last logical; TRUE if only extrapolation for the additional sample of size M
#' @param seed Set seed for randomness
#' 
#' @export
#' 
#' @details Draw samples from the posterior distribution of k_n + K^(n)_m, for different m,
#' for IBP with Gamma(a,b) mixture, and prior on a and b
extrapolation.GammaIBP <- function(object, M = 10, only_last = FALSE, seed = 1234) {
  
  # Extract the size of the observed sample "n" and the number of distinct features "Kn"
  feature_matrix <- object$feature_matrix
  n <- nrow(feature_matrix)
  Kn <- ncol(feature_matrix)
  
  # Perform the prediction
  extrapolation_list <- prediction_GammaIBP(object, n = n, Kn = Kn, M = M, only_last = only_last, seed = seed)
  
  return (extrapolation_list)
}
