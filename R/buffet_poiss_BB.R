##################################################################
# Poisson mixture of BB: sampling the features from the buffet generalization
##################################################################
# Input:
# alpha <0
# theta > -alpha
# n: number of individuals
# lambda > 0: parameter of Poisson mixing distribution

buffet_sampling_poiss_BB <- function(alpha, theta, n, lambda){
  
  
  n_dish <- rpois(1, -lambda*alpha/theta)
  
  features <- replicate(n_dish, list(c(1)))
  
  
  
  features[[1]] = c(features[[1]],2)
  
  
}