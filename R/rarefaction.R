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
#' (if n_reordering=1, then it is just the observed accumulation curve, with no reshuffling)
#' 
#' @param object  The \code{n x K}-dimensional binary matrix
#' @param n_reorderings number of reorderings of the data to average over
#' @param seed
#' 
#' @export
#' 
rarefaction.array <- function(object, n_reorderings = 1, seed = 1234) {
  
  feature_list <- convert_features_list(object)
  n <- nrow(object)
  
  if (n_reorderings == 1){
    
    rare_curve <- sapply(1:n, function(i) length(unique(unlist(feature_list[1:i]))) )
      
  } else {
    
    rare_curve <- matrix(NA, nrow = n_reorderings, ncol = n)
      
    for (j in 1:n_reorderings){
      
      f_list <- sample(feature_list)
      
      rare_curve[j, ] <- sapply(1:n, function(i) length(unique(unlist(f_list[1:i]))) )
    }
    
    rare_curve <- colMeans(rare_curve)
  }
  
  return(rare_curve)
  
}




#' Model-based rarefaction for the GibbsFA models
#'
#' @param object An object of class \code{GibbsFA}
#' @param seed Set seed for randomness
#' 
#' @export
#' 
#' @details Draw samples from the posterior distribution of K_n,
#' for BB with NB(n0,mu0) mixture
rarefaction.GibbsFA <- function(object, ...) {
  
  # Extract the size of the observed sample "n"
  feature_matrix <- object$feature_matrix
  n <- nrow(feature_matrix)

  # Perform the prediction
  if (class(object)[2] %in% c("PoissonBB", "NegBinBB", "GammaIBP")){
    rarefaction_list <- prediction(object, n = 0, Kn = 0, M = n, seed = seed)
    names(rarefaction_list) <- paste0("N = ", 1:n)
  } else {
    rarefaction_list <- prediction(object, n = 0, Kn = 0, M = n)
  }
  
  return (rarefaction_list)
}

