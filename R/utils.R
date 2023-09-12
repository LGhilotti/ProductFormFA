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
    #return (list("estimated" = est_new_features, "Umn" = 0))
    return (-1)
  }
  else {
    ratio <- abs(U_M_N - est_new_features)/U_M_N
  
    acc <- 1 - min(ratio, 1)
    
    #return (list("acc" = acc, "Umn" = U_M_N))
    return (acc)
  }
      
}



#' Metric of accuracy in prediction
#' 
#' 
#' @param obs_n [list] number of new observed
#' @param est_n [list] number of estimated new
#' @param obs_t [numeric] number of observed in training set
#' 
#' @export
#'
compute_accuracy <- function(obs_n, est_n, obs_t) {
  return (abs(obs_n - est_n)/obs_t)
}




#### Richness estimators in frequentist literature ############


#' Chao2 estimator
#' 
#' 
#' @param obs_n [list] number of new observed
#' @param est_n [list] number of estimated new
#' @param obs_t [numeric] number of observed in training set
#' 
#' @export
#'
chao2_estimator <- function(data_mat){
  
  # Compute total number of sites
  n <- nrow(data_mat)
  
  # Delete NA
  data_mat <- data_mat[, colSums(is.na(data_mat))==0]
  
  # Delete zero-columns
  data_mat <- data_mat[, colSums(data_mat)!=0]
  
  # Set K to be the observed number of features
  K <- ncol(data_mat)
  
  # Compute number of species contained in exactly k individuals
  counts <- colSums(data_mat)
  
  Q_1 <- sum(counts == 1)
  Q_2 <- sum(counts == 2)

  
  # Compute the statistic
  res <- K + (n-1)/n * (Q_1^2)/(2*Q_2)
  
  return (res)
  
}


#' Improved Chao2 estimator
#' 
#' 
#' @param obs_n [list] number of new observed
#' @param est_n [list] number of estimated new
#' @param obs_t [numeric] number of observed in training set
#' 
#' @export
#'
improved_chao2_estimator <- function(data_mat){
  
  # Compute total number of sites
  n <- nrow(data_mat)
  
  # Delete NA
  data_mat <- data_mat[, colSums(is.na(data_mat))==0]
  
  # Delete zero-columns
  data_mat <- data_mat[, colSums(data_mat)!=0]
  
  # Set K to be the observed number of features
  K <- ncol(data_mat)
  
  # Compute number of species contained in exactly k individuals
  counts <- colSums(data_mat)
  
  Q_1 <- sum(counts == 1)
  Q_2 <- sum(counts == 2)
  Q_3 <- sum(counts == 3)
  Q_4 <- sum(counts == 4)
  
  # Compute the statistic
  res <- K + (n-1)/n * (Q_1^2)/(2*Q_2) + (n-3)/(4*n) *Q_3/Q_4 * 
    max(0, Q_1 - (n-3)/(2*(n-1)) *Q_2*Q_3/ Q_4)
  
  return (res)
  
}


#' Beta-Binomial Chiu estimator
#' 
#' 
#' @param obs_n [list] number of new observed
#' @param est_n [list] number of estimated new
#' @param obs_t [numeric] number of observed in training set
#' 
#' @export
#'
beta_binomial_estimator <- function(data_mat){
  
  # Compute total number of sites
  n <- nrow(data_mat)
  
  # Delete NA
  data_mat <- data_mat[, colSums(is.na(data_mat))==0]
  
  # Delete zero-columns
  data_mat <- data_mat[, colSums(data_mat)!=0]
  
  # Set K to be the observed number of features
  K <- ncol(data_mat)
  
  # Compute number of species contained in exactly k individuals
  counts <- colSums(data_mat)
  
  Q_1 <- sum(counts == 1)
  Q_2 <- sum(counts == 2)
  Q_3 <- sum(counts == 3)
  
  # Compute Q_hat_0 
  if (Q_2 == 0){
    Q_hat_0 <- (n-1)/n * (Q_1*(Q_1 -1))/2
  } else {
    Q_hat_0 <- (n-1)/n * (Q_1^2)/(2*Q_2)
  }
  
  # Compute last term
  if (2*Q_2^2 / (3*Q_1*Q_3) <= 1){
    last <- 2 - max(0.5, 2*Q_2^2 / (3*Q_1*Q_3))
  } else {
    last <- 1
  }
  
  # Compute the statistic
  res <- K + Q_hat_0 * last
  
  return (res)
  
}