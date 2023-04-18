#' Credible intervals of Kmn for the BB with Poisson mixture
#'
#' This function computes the means and the credible intervals of Kmn for the BB with Poisson mixture
#'
#' @param alpha [numeric] value of alpha in product-form feature allocation
#' @param theta [numeric] value of theta in product-form feature allocation
#' @param m [numeric] dimension of the new sample to be observed
#' @param n [numeric] dimension of the already observed sample
#' @param lambda [numeric] Poisson hyperparameter 
#' @param lev [numeric] level of the credible intervals
#'
#' @export
#'
CI_Kmn_poiss_BB <- function(alpha, theta, m, n, lambda, lev) {
  
  means <- mean_kmn_all_poiss_BB(alpha, theta, m, n, lambda)
  ub <- qpois(lev + (1 - lev) / 2, means)
  lb <- qpois((1 - lev) / 2, means)

  return(list("means" = means,"ubs" = ub,"lbs" = lb))
}

#########################################################################

#' EB based on EFPF-max - BB with Poisson mixture
#'
#' This function returns the value of (alpha, theta, lambda) maximizing the 
#' EFPF for the given sample - BB with Poisson mixture
#'
#' @param n [integer] dimension of the observed sample
#' @param counts [numeric] vector of cardinalities for the observed features
#' @param pars_0 [numeric] Initialization for (alpha, theta, lambda) to 
#' be optimized
#' 
#' @export
EB_EFPF_poiss_BB <- function(n, counts, pars_0){
  
  # initialize s = theta + alpha
  pars_0[2] <- pars_0[2] + pars_0[1]
  
  # set constraints stricter so that function is limited
  res <- optim(par = pars_0, fn = neg_log_EFPF_poiss_BB_rep, n = n, counts = counts,
               method = "L-BFGS-B", lower = c(-Inf,0.1,0.1), upper = c(-0.1, Inf, Inf))
  
  sol <- res$par
  # convert to the theta parameter
  sol[2] <- sol[2] - sol[1]
  
  return (sol)
  
}


#########################################################################
######## Method of moments ##############################################
#########################################################################

#' EB based on method of moments - BB with Poisson mixture
#'
#' This function returns the value of (alpha, theta, lambda) obtained 
#' via the method of moments for the given sample - BB with Poisson mixture
#'
#' @param n [integer] dimension of the observed sample
#' @param ntrain [integer] dimension of training set 
#' @param val_rep [integer] number of samples in the validation averaged over for 
#' assessing the number of new features observed 
#' @param data_list [list] list of features by individual
#' @param pars_0 [numeric] Initialization for (alpha, theta, lambda) to 
#' be optimized
#' @param seed_id [numeric] Seed to be used
#' 
#' @export
EB_MM_poiss_BB <- function(n, ntrain, val_rep, data_list, pars_0, seed_id){
  
  set.seed(seed_id)
  # initialize s = theta + alpha
  pars_0[2] <- pars_0[2] + pars_0[1]
  
  # dimension of the validation set
  M <- n - ntrain
  
  # set constraints stricter so that function is limited
  res <- optim(par = pars_0, fn = mm_obj_poiss_BB_rep, ntest = M,
               ntrain=ntrain, val_rep = val_rep, data_list=data_list, seed_id = seed_id,
               method = "L-BFGS-B", lower = c(-Inf,0.1,0.1), upper = c(-0.1, Inf, Inf))
  
  sol <- res$par
  # convert to the theta parameter
  sol[2] <- sol[2] - sol[1]
  
  return (sol)
  
}


#########################################################################

#' MM objective function for BB with Poisson mixture with reparametrization
#'
#' @param pars [numeric] pars[1] = value of alpha in product-form feature allocation,
#' pars[2] = value of s = theta+alpha in product-form feature allocation,
#' pars[3] =  value of lambda - Poisson hyperparameter
#' 
#' @param ntest [integer] dimension of the validation set
#' @param ntrain [integer] dimension of the training set
#' @param val_rep [integer] number of samples in the validation averaged over for 
#' assessing the number of new features observed 
#' @param data_list [list] list of features by individual
#' @param seed_id [numeric] Seed to be used
#' 
#' @return value of the MM objective function to be minimized
#' 
#' @export
mm_obj_poiss_BB_rep <- function(pars, ntest, ntrain, val_rep, data_list, seed_id = seed_id){
  
  L <- ntrain + ntest
  alpha <- pars[1]
  s <- pars[2]
  lambda <- pars[3]
  theta <- s-alpha
  
  # training and test list
  train_list <- data_list[1:ntrain]
  test_list <- data_list[(ntrain+1):L]
  
  # features observed in the training set
  feat_train <- unique(unlist(train_list))
  
  # storing the sum
  sum <- 0
  
  for (m in 1:ntest){
    # features observed in the sampled validation
    nfeat_new_obs_avg <- 0
    set.seed(seed_id)
    for (r in 1:val_rep){
      idx_val_rep <- sample(1:ntest, size = m) 
      feat_val <- unique(unlist( test_list[idx_val_rep] ))
      
      nfeat_new_obs <- length(setdiff(feat_val, feat_train))
      
      nfeat_new_obs_avg <- nfeat_new_obs_avg + nfeat_new_obs
    }
    # average number of observed new features (Umn observed)
    nfeat_new_obs_avg <- nfeat_new_obs_avg/val_rep
    
    # predicted number of new features
    nfeat_new_pred <- mean_kmn_poiss_BB(alpha, theta, m, ntrain, lambda)
    
    sum <- sum + (nfeat_new_pred - nfeat_new_obs_avg)**2
  }
  
  return (sum)
    
}


#######################################################################
############ Method of moments - cross-validation #####################
#######################################################################

#########################################################################

#' EB based on method of moments with cross-validation - BB with Poisson mixture
#'
#' This function returns the value of (alpha, theta, lambda) obtained 
#' via the method of moments with cross-validation for the given sample -
#'  BB with Poisson mixture
#'
#' @param n [integer] dimension of the observed sample
#' @param nfold [integer] number of folds 
#' @param obs_sample [list] list with the whole sample of features, where
#' $features contains the simulated features for each customer,
#' $num_new contains the number of new features selected for each customer
#' $counts contains the counts for the observed features
#' @param pars_0 [numeric] Initialization for (alpha, theta, lambda) to 
#' be optimized
#' @param seed_id [numeric] Seed to be used
#' 
#' @export
EB_MM_cv_poiss_BB <- function(n, nfold, obs_sample, pars_0, seed_id){
  
  set.seed(seed_id)
  # initialize s = theta + alpha
  pars_0[2] <- pars_0[2] + pars_0[1]
  
  # reord_idx is the new order of the data to consider
  reord_idx <- sample(1:n)
  # card_folds is a nfold-dimensional vector containing the dimensions of the folds
  card_folds <- rep(floor(n/nfold), nfold)
  rem <- n - nfold*floor(n/nfold)
  if (rem > 0) card_folds[1:rem] = card_folds[1:rem] +1
  
  # set constraints stricter so that function is limited
  res <- optim(par = pars_0, fn = mm_obj_cv_poiss_BB_rep, reord_idx = reord_idx,
               card_folds = card_folds, obs_sample = obs_sample, 
               method = "L-BFGS-B", lower = c(-Inf,0.1,0.1), upper = c(-0.1, Inf, Inf))
  
  sol <- res$par
  # convert to the theta parameter
  sol[2] <- sol[2] - sol[1]
  
  return (sol)
  
}



#########################################################################

#' MM objective function with cross-validation for BB with Poisson mixture 
#' with reparametrization
#'
#' @param pars [numeric] pars[1] = value of alpha in product-form feature allocation,
#' pars[2] = value of s = theta+alpha in product-form feature allocation,
#' pars[3] =  value of lambda - Poisson hyperparameter
#' 
#' @param reord_idx [numeric] order of the observations to be considered
#' @param card_folds [numeric] cardinalities of the folds
#' @param obs_sample [list] list with the whole sample of features, where
#' $features contains the simulated features for each customer,
#' $num_new contains the number of new features selected for each customer
#' $counts contains the counts for the observed features#' 
#' 
#' @return value of the MM objective function with cross-validation to be
#'  minimized
#' 
#' @export
mm_obj_cv_poiss_BB_rep <- function(pars, reord_idx, card_folds, obs_sample){
  
  nfold <- length(card_folds)
  
  fun <- 0
  
  # store the observation indexes and the observed features for each fold
  idx_fold <- vector("list",nfold)
  feat_by_fold <- vector("list",nfold)
  
  idx_fold[[1]] <- reord_idx[1:card_folds[1]]
  feat_by_fold[[1]] <- unique(unlist( obs_sample$features[idx_fold[[1]]] ))
  for (j in 2:nfold){
    # store the features observed in each fold
    idx_fold[[j]] <- reord_idx[(sum(card_folds[1:(j-1)])+1) : sum(card_folds[1:j])]
    feat_by_fold[[j]] <- unique(unlist( obs_sample$features[idx_fold[[j]]] ))
  }
  
  # perform the single evaluation of the cross-validation
  for (j in 1:nfold){
    # features observed in the folds different from j
    feat_others <- unique(unlist(feat_by_fold[-j]))
    # number of features contained in fold j and not in the others
    num_new_fold_m <- rep(0, card_folds[j])
    # features observed in individuals in fold, ordered as reord_idx
    feat_in_fold <- obs_sample$features[idx_fold[[j]]]
    for (i in 1:card_folds[j]){
      feat_in_fold_till_m <- unique(unlist(feat_in_fold[1:i]))
      num_new_fold_m <- setdiff()
    }
    num_new_fold <- length(setdiff(feat_by_fold[[j]], feat_others))
    
  }
  
}
