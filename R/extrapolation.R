#
# Functions for computing the extrapolation curve, i.e. the statistic k_n + K^(n)_m, a posteriori.
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



#' Model-based extrapolation for GibbsFA models
#'
#' @param object An object of class \code{GibbsFA}
#' @param M Maximum size of future sample to extrapolate
#' @param only_last logical; TRUE if only extrapolation for the additional sample of size M
#' @param seed Set seed for randomness
#' 
#' @export
#' 
#' @details Draw samples from the posterior distribution of k_n + K^(n)_m, for different m,
#' for the GibbsFA models
extrapolation.GibbsFA <- function(object, M = 10, only_last = FALSE, ...) {
  
  # Extract the size of the observed sample "n" and the number of distinct features "Kn"
  feature_matrix <- object$feature_matrix
  n <- nrow(feature_matrix)
  Kn <- ncol(feature_matrix)

  # Perform the prediction
  if (class(object)[2] %in% c("PoissonBB", "NegBinBB", "GammaIBP")){
    extrapolation_list <- prediction(object, n = n, Kn = Kn, M = M, only_last = only_last, seed = seed)
  } else{
    extrapolation_list <- prediction(object, n = n, Kn = Kn, M = M, only_last = only_last)
  }
  
  return (extrapolation_list)
}

