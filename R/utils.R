#' Create features matrix
#' 
#' This function allows to create the matrix of the 
#' observed features in order-of-appearance, with 1 if present, 0 if not
#'
#' @param feat_list [list] list containing the observed features
#'
#' @export
#'
create_features_matrix <- function(feat_list){
  
  n_individuals <- length(feat_list)
  mat <- matrix(0,n_individuals,max(unlist(feat_list)))
  
  for (i in 1:n_individuals){
    mat[i,feat_list[[i]]] <- 1 
  }
  
  return(mat)
  
}



#' Convert matrix of 0/1 to list of individuals with observed features
#' 
#' This function allows to create the list of individuals with observed features 
#' from the matrix of 0/1
#' 
#' @param mat [integer] matrix of 0/1 to be converted
#'
#' @export
#'
create_features_list <- function(mat){
  
  feat_list <- vector("list", nrow(mat))
  
  for (i in 1:nrow(mat)){
    feat_list[[i]] <- which(mat[i,]==1, arr.ind = TRUE)
  }
  
  return (feat_list)
}


#' Percentage accuracy in prediction
#' 
#' This function computes the percentage accuracy of the proposed method on
#' the test set
#' 
#' @param train_list [list] list containing the training observations
#' @param test_list [list] list containing the test observations
#' @param est_new_features [numeric] estimated number of new features in the test set
#' (computed with N=#train_list and M=#test_list)
#'
#' @return the percentage accuracy of the estimate on the test set
#' 
#' @export
#'
perc_accuracy <- function(train_list, test_list, est_new_features){
  
  feat_train <- unique(unlist(train_list))
  feat_test <- unique(unlist(test_list))
  obs_new_features <- setdiff(feat_test, feat_train)
  U_M_N <- length(obs_new_features)
  
  if (U_M_N==0){
    return (list("estimated" = est_new_features, "Umn" = 0))
  }
  else {
    ratio <- abs(U_M_N - est_new_features)/U_M_N
  
    acc <- 1 - min(ratio, 1)
    
    return (list("acc" = acc, "Umn" = U_M_N))
  }
      
}