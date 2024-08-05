#
# Functions for computing the K_n_r curves.
# Available for both mixtures of Beta-Bernoulli and Gamma IBP.
#

#' Model-based K_n_r
#' @param object An object of class \code{GibbsFA}
#' @param ... Additional parameters
#' 
#' @export
K_n_r <- function(object, ...) {
  
  UseMethod("K_n_r", object)
}


#' Sample-based K_n_r curves for a vector of species abundances 
#' (if n_reordering=1, then it is just the observed curves, with no reshuffling)
#' 
#' @param object  The \code{n x K}-dimensional binary matrix
#' @param n_reorderings number of reorderings of the data to average over
#' @param seed
#' 
#' @export
#' 
K_n_r.array <- function(object, n_reorderings = 1, seed = 1234) {
  
  feature_list <- convert_features_list(object)
  n <- nrow(object)
  
  K_n_r_curves <- vector("list", n)
  # element i is a i-dim vector of K_i_r, with r =1,...,i
  
  if (n_reorderings == 1){
    
   curves <- lapply(1:n, function(i) table(table(unlist(feature_list[1:i]))) )
   
   K_n_r_curves <- lapply(1:n, function(i) sapply(1:i, function(r){
     if (r %in% names(curves[[i]])){
       return(unname(curves[[i]][names(curves[[i]]) == r]))
     } else {
       return(0)
     }
     }))
     

  } else {
    
    K_n_r_curves <- lapply(1:n, function(i) rep(0, i) )
    
    for (j in 1:n_reorderings){
      
      f_list <- sample(feature_list)
      
      curves <- lapply(1:n, function(i) table(table(unlist(f_list[1:i]))) )
      
      K_n_r_curves <- mapply(FUN = function(x,y) x + y, K_n_r_curves,
                             lapply(1:n, function(i) sapply(1:i, function(r){
                               if (r %in% names(curves[[i]])){
                                 return(unname(curves[[i]][names(curves[[i]]) == r]))
                               } else {
                                 return(0)
                               }
                             })))
      
    }
    
    K_n_r_curves <- lapply(K_n_r_curves, function(list_i) list_i/ n_reorderings)
  
  }
  
  names(K_n_r_curves) <- paste0("N = ", 1:n)
  
  return(K_n_r_curves)
  
}




#' Model-based K_n_r for the GibbsFA models
#'
#' @param object An object of class \code{GibbsFA}
#' @param seed Set seed for randomness
#' 
#' @export
#' 
#' @details Perform estimate on K_n_r
K_n_r.GibbsFA <- function(object, ...) {
  
  # Extract the size of the observed sample "n"
  feature_matrix <- object$feature_matrix
  n <- nrow(feature_matrix)
  
  # Perform the estimation
  if (class(object)[2] %in% c("PoissonBB", "NegBinBB", "GammaIBP")){
    K_n_r_list <- prediction_K_n_r(object, n, seed = seed)
    names(K_n_r_list) <- paste0("N = ", 1:n)
  } else {
    K_n_r_list <- prediction_K_n_r(object, n)
    names(K_n_r_list) <- paste0("N = ", 1:n)
  }
  
  return (K_n_r_list)
}



#' Model-based extrapolation
#' @param object An object of class \code{GibbsFA}
#' @param ... Additional parameters
#' 
#' @export
prediction_K_n_r <- function(object, ...) {
  
  UseMethod("prediction_K_n_r", object)
}


#' Model-based K_n_r for BB with Poisson(lambda) mixture - EB version
#'
#' @param object An object of class \code{GibbsFA, PoissonBB_eb}
#' @param n Size of the hypothetical observed sample
#' 
#' @details Return list of lambdas sequence of the K_n_r, for different n,
#' for BB with Poisson(lambda) mixture
prediction_K_n_r.PoissonBB_eb <- function(object, n) {
  
  # Extract the parameters from object
  alpha <- object$alpha
  theta <- object$theta
  lambda <- object$lambda
  
  poiss_par_list <- lapply(1:n, function(i) 
    list("lambda_est" = lambda*(-alpha)*choose(i, 1:i)*exp(
      lgamma(1:i - alpha) - lgamma(1 - alpha) +
        lgamma(theta + alpha + i - 1:i) - lgamma(theta + alpha) -
        lgamma(theta + i) + lgamma(theta)))
    ) 
  
  return (poiss_par_list)
}



#' Model-based K_n_r for BB with NB(n0,mu0) mixture - EB version
#'
#' @param object An object of class \code{GibbsFA, NegBinBB}
#' @param n Size of the hypothetical observed sample
#' 
#' @details Return list of the (n0,mu0)s sequences of the NB of K_n_r, for different n,
#' for BB with NB(n0,mu0) mixture
prediction_K_n_r.NegBinBB_eb <- function(object, n) {
  
  # Extract the parameters from object
  alpha <- object$alpha
  theta <- object$theta
  n0 <- object$n0
  mu0 <- object$mu0
  
  nb_pars_list <- lapply(1:n, function(i) 
    list("mu0_est" = mu0*(-alpha)*choose(i, 1:i)*exp(
      lgamma(1:i - alpha) - lgamma(1 - alpha) +
        lgamma(theta + alpha + i - 1:i) - lgamma(theta + alpha) -
        lgamma(theta + i) + lgamma(theta)),
      "n0_est" = rep(n0, i))
  ) 
  
  return (nb_pars_list)
}


#' Model-based K_n_r for IBP with Gamma(a,b) mixture - EB version
#'
#' @param object An object of class \code{GibbsFA, GammaIBP}
#' @param n Size of the hypothetical observed sample
#' 
#' @details Return list of the (n0,mu0) sequences of the NB for K_n_r, 
#' for different n, for IBP with Gamma(a,b) mixture
prediction_K_n_r.GammaIBP_eb <- function(object, n) {
  
  # Extract the parameters from object
  alpha <- object$alpha
  theta <- object$theta
  a <- object$a
  b <- object$b
  
  nb_pars_list <- lapply(1:n, function(i) 
    list("mu0_est" = a/b *choose(i, 1:i)*exp(
      lgamma(theta + 1) - lgamma(theta + n) +
      lgamma(1:i - alpha) - lgamma(1 - alpha) +
        lgamma(theta + alpha + i - 1:i) - lgamma(theta + alpha)),
      "n0_est" = rep(a, i))
  ) 
  
  return (nb_pars_list)
  
}



prediction_K_n_r.GammaIBP <- function(object, n, seed) {
  stop("Not implemented")
}

prediction_K_n_r.PoissonBB <- function(object, n, seed) {
  stop("Not implemented")
}

prediction_K_n_r.NegBinBB <- function(object, n, seed) {
  stop("Not implemented")
}


ev_K_n_r_GammaIBP <- function(alpha, theta, a, b, n){
  
  return( a/b *choose(n, 1:n)*exp(
    lgamma(theta + 1) - lgamma(theta + n) +
      lgamma(1:n - alpha) - lgamma(1 - alpha) +
      lgamma(theta + alpha + n - 1:n) - lgamma(theta + alpha))
    )
  
}

ev_K_n_r_PoissonBB <- function(alpha, theta, lambda, n){
  
  return( lambda*(-alpha)*choose(n, 1:n)*exp(
    lgamma(1:n - alpha) - lgamma(1 - alpha) +
      lgamma(theta + alpha + n - 1:n) - lgamma(theta + alpha) -
      lgamma(theta + n) + lgamma(theta))
  )
  
}

ev_K_n_r_NegBinBB <- function(alpha, theta, n0, mu0, n){
  
  return( mu0*(-alpha)*choose(n, 1:n)*exp(
    lgamma(1:n - alpha) - lgamma(1 - alpha) +
      lgamma(theta + alpha + n - 1:n) - lgamma(theta + alpha) -
      lgamma(theta + n) + lgamma(theta))
  )
  
}