#' Create features matrix
#' 
#' This function allows to create the matrix of the 
#' observed features in order-of-appearance, with 1 if present, 0 if not
#'
#' @param feature_list [list] list containing the observed features
#'
#' @export
#'
convert_features_matrix <- function(feature_list){
  
  n_individuals <- length(feature_list)
  mat <- matrix(0,n_individuals,max(unlist(feature_list)))
  
  for (i in 1:n_individuals){
    mat[i,feature_list[[i]]] <- 1 
  }
  
  return(mat)
  
}



#' Convert matrix of 0/1 to list of individuals with observed features
#' 
#' This function allows to create the list of individuals with observed features 
#' from the matrix of 0/1
#' 
#' @param feature_matrix [integer] matrix of 0/1 to be converted
#'
#' @export
#'
convert_features_list <- function(feature_matrix){
  
  feat_list <- vector("list", nrow(feature_matrix))
  
  for (i in 1:nrow(feature_matrix)){
    feat_list[[i]] <- which(feature_matrix[i,]==1, arr.ind = TRUE)
  }
  
  return (feat_list)
}
