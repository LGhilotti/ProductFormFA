##################################################################
###### Simulation of features from Zipf ##########################
##################################################################

#' Simulate observations with features from Zipf distribution
#'
#' This function draws samples of features/observations from Zipf distribution
#'
#' @param n [integer] number of observations to extract
#' @param K [integer] maximum number of observable features
#' @param xi [numeric] probability factor
#'
#' @export
#'
rzipf <- function(n, K, xi, seed = 1234){
  
  set.seed(seed)
  
  prob_vec <- sapply(1:K, function(x) (1/(x+1))**xi)
  
  zipf_sample <- matrix(, nrow = n, ncol = K)
  
  for (i in 1:n){
    zipf_sample[i,] <- rbinom(n = K, size = 1, prob = prob_vec)
  }
  
  return (zipf_sample)
    
}