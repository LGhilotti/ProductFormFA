#' Create features matrix
#' 
#' This function allows to create the matrix of the 
#' observed features in order-of-appearance, with 1 if present, 0 if not
#'
#' @param feat_list [list] feat_list$features contains the observed features, 
#' feat_list$num_feat contains the total number of features
#' 
#'
#' @export
#'
create_features_matrix <- function(feat_list){
  
  n_individuals <- length(feat_list$features)
  mat <- matrix(0,n_individuals,feat_list$num_feat)
  
  for (i in 1:n_individuals){
    mat[i,feat_list$features[[i]]] <- 1 
  }
  
  return(mat)
  
}