
rm(list=ls())
library(ProductFormFA)
library(tidyverse)
library(scales)
library(ggpubr)
library(ggthemes)

source("R_script_paper/Routine_Chao.R")
source("R_script_paper/utils.R")

# Functions for fitting the SB-SP model -----

#################### Functions for SB-SP model #################
compute_der_log_full_sb_sp <- function(alpha_hat, 
                                       c, beta, Gam, n, K, counts,
                                       a_alpha, b_alpha){
  
  ealpha <- exp(alpha_hat)
  
  p2 <- exp( lgamma((1:n)) - 
               lgamma(1 - ealpha/(1+ealpha) + (1:n)) + lgamma(2 - ealpha/(1+ealpha) ) )
  
  # contribution of the Gamma-IBP
  dalpha_hat <- a_alpha + ealpha/((1 + ealpha)^2) *
    ( K*(digamma(1 - ealpha/(1+ealpha) + n ) - digamma(2 - ealpha/(1+ealpha) ) ) -
        Gam*sum( p2 * (digamma(1 - ealpha/(1+ealpha) + (1:n)) - digamma(2 - ealpha/(1+ealpha) ))) -
        sum(digamma(1/(1 + ealpha) + counts - 1)) + 
        K*digamma(1/(1+ealpha)) - 
        (a_alpha + b_alpha)*(1+ ealpha))
  
  # due to the fact alpha is in the prior on Gam
  dalpha_hat <- dalpha_hat - (c + 1) + Gam*beta/ealpha
  
  return (dalpha_hat)
  
}


compute_log_ratio_full_sb_sp <- function(alpha_hat_prop,
                                         alpha_hat_curr, 
                                         c, beta, Gam, n, K, counts,
                                         a_alpha, b_alpha){
  
  ealpha_prop <- exp(alpha_hat_prop)
  ealpha_curr <- exp(alpha_hat_curr)
  
  p2_prop <- exp( lgamma( (1:n) ) - 
                    lgamma(1 - ealpha_prop/(1+ealpha_prop) + (1:n)) + 
                    lgamma(2 - ealpha_prop/(1+ealpha_prop) ) )
  
  p2_curr <- exp( lgamma( (1:n) ) - 
                    lgamma(1 - ealpha_curr/(1+ealpha_curr) + (1:n)) + 
                    lgamma(2 - ealpha_curr/(1+ealpha_curr) ) )
  
  v <- lgamma(1/(1 + ealpha_prop) + counts - 1) - lgamma(1/(1 + ealpha_prop)) -
    lgamma(1/(1 + ealpha_curr) + counts - 1) + lgamma(1/(1 + ealpha_curr)) 
  
  # log-ratio of the Gamma-IBP contribution
  res <- K*(lgamma(2 - ealpha_prop/(1+ealpha_prop) ) - 
              lgamma(1 - ealpha_prop/(1+ealpha_prop) + n) -
              lgamma(2 - ealpha_curr/(1+ealpha_curr) ) + 
              lgamma(1 - ealpha_curr/(1+ealpha_curr) + n)) -
    Gam*sum(p2_prop - p2_curr) +
    sum(v) + 
    a_alpha*(alpha_hat_prop - alpha_hat_curr) - 
    (a_alpha + b_alpha)*(log(1 + ealpha_prop) - log(1 + ealpha_curr))
  
  # due to the fact alpha is in prior on Gam
  res <- res + (c+1)*(alpha_hat_curr - alpha_hat_prop) + 
    Gam*beta*(1/ealpha_curr - 1/ealpha_prop)
  
  
  return (res)
  
}


compute_log_ratio_q_sb_sp <- function(alpha_hat_prop,
                                      alpha_hat_curr, tau,
                                      der_log_full_curr, der_log_full_prop){
  
  norm_prop <- (alpha_hat_curr - alpha_hat_prop - tau*der_log_full_prop)^2
  
  norm_curr <- (alpha_hat_prop - alpha_hat_curr - tau*der_log_full_curr)^2
  
  res <- 1/(4*tau)*(norm_curr - norm_prop)
  
  return (res)
  
}


sampler_SB_SP <- function(Z,
                          c_0, beta_0, alpha_0,
                          p, r, t, a_alpha, b_alpha, 
                          tau, S, n_burnin, thin, seed){
  
  set.seed(seed)
  
  # Compute total number of sites
  n <- nrow(Z)
  
  # Delete NA
  Z <- Z[, colSums(is.na(Z))==0]
  
  # Delete zero-columns
  Z <- Z[, colSums(Z)!=0]
  
  # Set K to be the observed number of features
  K <- ncol(Z)
  
  # Compute vector of counts
  counts <- colSums(Z)
  
  ############## Gibbs-sampler ##########################
  
  # Define structure to store parameters along the iterations
  number_saved_iterations <- (S - n_burnin)/thin + 1
  c_vec <- vector(length = number_saved_iterations)
  beta_vec <- vector(length = number_saved_iterations)
  alpha_vec <- vector(length = number_saved_iterations)
  
  # Set initial values
  c <- c_0
  beta <- beta_0
  alpha <- alpha_0
  
  # index saved iterations (after burn-in and thinning satisfied)
  l <- 1
  
  for (q in 1:S){
    
    ################################################################
    ################# Draw c, beta, alpha | Z   ###################
    ###############################################################
    
    # Update Gam | Z, alpha, c, beta
    
    gamma_a_1ma_n <- sum(exp(lgamma((1:n)) - 
                               lgamma(1 - alpha + (1:n)) + lgamma(2 - alpha) ) )
    
    Gam <- rgamma(1, shape = K + c + 1, rate = beta*(1-alpha)/alpha + gamma_a_1ma_n)
    
    # Update c | Gam, Z, alpha, beta
    
    c <- rpois(1, beta*Gam*(1-p)*(1-alpha)/alpha )
      
    # Update beta | Gam, Z, alpha, c
    
    beta <- rgamma(1, shape = c + 1 + r, rate = Gam*(1-alpha)/alpha + t )
  
    
    ################################################################
    ############# Draw alpha | Z, Gam, c, beta ##################
    ###############################################################
    
    # In order to update alpha, we update alpha_hat, 
    # defined as transformation of alpha
    
    ### Current values for alpha_hat
    alpha_hat_curr <- log(alpha/(1-alpha))
    
    ### Propose values for alpha_hat
    # Compute the gradient of log-full conditional for MALA 
    der_log_full_curr <- compute_der_log_full_sb_sp(alpha_hat_curr, 
                                                    c, beta, Gam, n, K, counts,
                                                    a_alpha, b_alpha)
    
    alpha_hat_prop <- alpha_hat_curr + tau*der_log_full_curr + sqrt(2*tau)*rnorm(1)
    
    
    ### Acceptance probability 
    # Compute the log ratio of the full-cond in prop point and curr point
    log_ratio_full <- compute_log_ratio_full_sb_sp(alpha_hat_prop,
                                                   alpha_hat_curr, 
                                                   c, beta, Gam, n, K, counts,
                                                   a_alpha, b_alpha)
    
    # Compute the log ratio of the terms related to the proposal q
    der_log_full_prop <- compute_der_log_full_sb_sp(alpha_hat_prop,
                                                    c, beta, Gam, n, K, counts,
                                                    a_alpha, b_alpha)
    
    log_ratio_q <- compute_log_ratio_q_sb_sp(alpha_hat_prop,
                                             alpha_hat_curr, tau,
                                             der_log_full_curr, der_log_full_prop)
    
    # Compute acceptance probability
    log_acc_prob <- log_ratio_full + log_ratio_q
    acc_prob <- min(1, exp(log_acc_prob))
    
    # Decide if accept or not the new alpha
    if (runif(1) < acc_prob){ # accept
      alpha <- exp(alpha_hat_prop)/(1 + exp(alpha_hat_prop))
    }
    
    
    #######################################################################
    
    # Store parameters if burn-in is over and once every "thin" iteration
    if ((q > n_burnin) & (q %% thin == 0) ){
      print(paste0("iteration: ", q))
      
      c_vec[l] <- c
      beta_vec[l] <- beta
      alpha_vec[l] <- alpha
      
      l <- l+1
    }
    
  }
  
  c_vec <- c_vec[1:(l-1)]
  beta_vec <- beta_vec[1:(l-1)]
  alpha_vec <- alpha_vec[1:(l-1)]
  
  return (list("c_chain" = c_vec, "beta_chain" = beta_vec, 
               "alpha_chain" = alpha_vec))
  
}



# Function to generate binary matrix according to the polynomial mechanism ----

generate_data <- function(xi, n, H, seed = 1234){
  
  set.seed(seed)
  
  pis <- 1/((1:H)^xi)
  
  data_mat <- matrix(rbinom(n*H, size = 1, prob = rep(pis, n)),
                     nrow = n, ncol = H, byrow = T )
  
  return(data_mat)
  
}


# Single dataset: function to produce fit and estimate on specific polynomial scenario ------

fit_estimate_polynomial_scenario_singledataset <- function(xi,
                                                           init_obj_PoissonBB, mcmcparams_obj_PoissonBB, prior_obj_PoissonBB,
                                                           init_obj_NegBinBB, mcmcparams_obj_NegBinBB, prior_obj_NegBinBB, c_fr,
                                                           init_obj_GammaIBP, mcmcparams_obj_GammaIBP, prior_obj_GammaIBP,
                                                           seed = 1234){
  
  
  # Set training dimensions and dimension of the whole sample 
  Ns <- c(20, 40, 80)
  L <- 500
  # Set maximum number of features
  H <- 10^6
  
  # Generate data 
  data_mat <- generate_data(xi = xi, n = L, H = H, seed = seed)
  
  # Structures to save fit and estimates 
  
  # List to store the MCMC chains for PoissonBB, NegBinBB and GammaIBP, under different settings 
  list_params_PoissonBB <- vector(mode = "list")
  list_params_NegBinBB <-  vector(mode = "list")
  list_params_GammaIBP <-  vector(mode = "list")
  
  # List to store extrapolation draws for PoissonBB, NegBinBB and GammaIBP, for different settings
  list_extr_PoissonBB <- vector(mode = "list")
  list_extr_NegBinBB <-  vector(mode = "list")
  list_extr_GammaIBP <-  vector(mode = "list")
  
  # List of Nbar_emp for the different training sets
  list_Nbar_emp <- vector(mode = "list")
  
  # List to store SB-SP estimates
  list_extr_SBSP <- vector(mode="list")
  # List to store smoothed GT's estimates
  list_extr_GT <- vector(mode="list")
  
  
  
  # Loop over the different training set dimensions 
  
  for (j in 1:length(Ns)){
    
    n_train <- Ns[j]
    M <- L - n_train
    
    train_mat <- data_mat[1:n_train,]
    train_mat <- train_mat[, colSums(is.na(train_mat))==0]
    train_mat <- train_mat[, colSums(train_mat)!=0]
    
    lab_comb_bb <- paste0("n_train.",n_train,":Nbar.emp")
    lab_comb_ibp <- paste("n_train", n_train, sep = ".")
    
    # Empirical estimate of E(N) is obtained by Chiu
    Nbar_emp <- beta_binomial_estimator(train_mat)
    list_Nbar_emp[[paste0("n_train.",n_train)]] <- Nbar_emp
    
    # Fit the models
    # PoissonBB
    prior_obj_PoissonBB$lambda <- Nbar_emp
    PoissonBB_fit <- GibbsFA(feature_matrix = train_mat, 
                             model = "PoissonBB", 
                             prior = prior_obj_PoissonBB, 
                             initialization = init_obj_PoissonBB, 
                             mcmcparams = mcmcparams_obj_PoissonBB)
    # NegBinBB
    prior_obj_NegBinBB$n0 <- Nbar_emp/(c_fr - 1)
    prior_obj_NegBinBB$mu0 <- 1/c_fr
    NegBinBB_fit <- GibbsFA(feature_matrix = train_mat, 
                            model = "NegBinBB", 
                            prior = prior_obj_NegBinBB, 
                            initialization = init_obj_NegBinBB, 
                            mcmcparams = mcmcparams_obj_NegBinBB)    
    # GammaIBP
    GammaIBP_fit <- GibbsFA(feature_matrix = train_mat, 
                            model = "GammaIBP", 
                            prior = prior_obj_GammaIBP, 
                            initialization = init_obj_GammaIBP, 
                            mcmcparams = mcmcparams_obj_GammaIBP)  
    
    
    # Fill the structures of MCMC chains of the parameters
    list_params_PoissonBB[[lab_comb_bb]] <- list("alpha" = PoissonBB_fit$alpha_chain, 
                                                 "theta" = PoissonBB_fit$theta_chain)
    
    list_params_NegBinBB[[lab_comb_bb]] <- list("alpha" = NegBinBB_fit$alpha_chain, 
                                                "theta" = NegBinBB_fit$theta_chain)
    
    list_params_GammaIBP[[lab_comb_ibp]] <- list("alpha" = GammaIBP_fit$alpha_chain, 
                                                 "theta" = GammaIBP_fit$theta_chain,
                                                 "a" = GammaIBP_fit$a_chain,
                                                 "b" = GammaIBP_fit$b_chain)
    
    
    # Fill the extrapolation structures
    list_extr_PoissonBB[[lab_comb_bb]] <- extrapolation(object = PoissonBB_fit, M = M)
    list_extr_NegBinBB[[lab_comb_bb]] <- extrapolation(object = NegBinBB_fit, M = M)
    list_extr_GammaIBP[[lab_comb_ibp]] <- extrapolation(object = GammaIBP_fit, M = M)
    
    
    # Competitors (actually, SB-SP is special case of GammaIBP)
    
    # A) SB-SP
    
    output_sp <- sampler_SB_SP(Z = train_mat,
                              c_0 = 10, beta_0 = 10, alpha_0 = 0.5,
                              p = 0.05, r = 0.1, t = 0.01, a_alpha = 2, b_alpha = 2,
                              tau = 0.005, S = 100, n_burnin = 10, thin = 2, seed = seed)
    
    SBSP_fit <- list("feature_matrix" = train_mat,
                     "alpha_chain" = output_sp$alpha_chain,
                     "theta_chain" = 1 - output_sp$alpha_chain,
                     "a_chain" = output_sp$c_chain +1,
                     "b_chain" = output_sp$beta_chain*(1-output_sp$alpha_chain)/output_sp$alpha_chain)
    
    class(SBSP_fit) <- c("GibbsFA", "GammaIBP")
    
    list_extr_SBSP[[lab_comb_ibp]] <- extrapolation(object = SBSP_fit, M = M)
    
    
    # B) Smoothed Good-Toulmin extrapolation
    
    # Compute SFS vector and CTS vector
    sfs <- tabulate(colSums(train_mat))
    cts <- sapply(2:n_train, function(n) ncol(train_mat[1:n,colSums(train_mat[1:n,]) > 0])   )
    cts <- c(0, sum(train_mat[1,]) , cts)
    
    list_extr_GT[[lab_comb_ibp]] <- predict_good_toulmin(n_train, M, sfs, cts, alternative = 0)$preds
    
  }
  
  # Save the entire workspace related to the mechanism just performed
  save(list = ls(all.names = TRUE), file =  paste0("R_script_paper/poly_",xi,"_fit_estimate_singledataset.RData"))
}


# EB Single dataset: function to produce fit and estimate on specific ecological scenario ------

eb_fit_estimate_polynomial_scenario_singledataset <- function(xi, 
                                                              var_fct,
                                                              eb_init_GammaIBP, eb_known_GammaIBP,
                                                              seed = 1234){
  
  if (xi <= 0){
    stop("Invalid generating mechanism.")
  }
  
  # Set training dimensions and dimension of the whole sample 
  Ns <- c(30, 60, 120) #c(20, 40, 80) 
  L <- 600
  # Set maximum number of features
  H <- 10^6
  
  # Generate data 
  data_mat <- generate_data(xi = xi, n = L, H = H, seed = seed)
  
  # Structures to save fit and estimates 
  
  # List to store the fitted models for PoissonBB, NegBinBB and GammaIBP, under different settings 
  list_eb_MM_cens_fit_PoissonBB <- vector(mode = "list")
  list_eb_MM_cens_fit_NegBinBB <-  vector(mode = "list")
  list_eb_EFPF_fit_GammaIBP <-  vector(mode = "list")
  list_eb_MM_fit_GammaIBP <-  vector(mode = "list")
  
  # Loop over the different training set dimensions 
  
  for (j in 1:length(Ns)){
    
    n_train <- Ns[j]
    M <- L - n_train
    
    train_mat <- data_mat[1:n_train,]
    
    lab_comb_bb <- paste0("n_train.",n_train,":Nbar.emp")
    lab_comb_ibp <- paste("n_train", n_train, sep = ".")
    
    # Fit the models
    # # PoissonBB
    # eb_MM_cens_PoissonBB_fit <- GibbsFA_eb(feature_matrix = train_mat, 
    #                                        model = "PoissonBB", 
    #                                        type = "MM_censored")
    # 
    # # NegBinBB
    # eb_MM_cens_NegBinBB_fit <- GibbsFA_eb(feature_matrix = train_mat, 
    #                                       model = "NegBinBB", 
    #                                       type = "MM_censored", var_fct = var_fct)
    # 
    # GammaIBP
    eb_params_obj_GammaIBP <- eb_params(model = "GammaIBP",
                                        init = eb_init_GammaIBP,
                                        known = eb_known_GammaIBP)
    
    eb_EFPF_GammaIBP_fit <- GibbsFA_eb(feature_matrix = train_mat,
                                     model = "GammaIBP", type = "EFPF",
                                     eb_params =  eb_params_obj_GammaIBP)
    
    eb_MM_GammaIBP_fit <- GibbsFA_eb(feature_matrix = train_mat,
                                     model = "GammaIBP", type = "MM",
                                     eb_params =  eb_params_obj_GammaIBP)
    
    
    # Fill the structures of MCMC chains of the parameters
    list_eb_MM_cens_fit_PoissonBB[[lab_comb_bb]] <- eb_MM_cens_PoissonBB_fit
    
    list_eb_MM_cens_fit_NegBinBB[[lab_comb_bb]] <- eb_MM_cens_NegBinBB_fit
    
    list_eb_EFPF_fit_GammaIBP[[lab_comb_bb]] <- eb_EFPF_GammaIBP_fit
    list_eb_MM_fit_GammaIBP[[lab_comb_bb]] <- eb_MM_GammaIBP_fit
    
    
    # # Competitors
    # # B) Smoothed Good-Toulmin extrapolation
    # 
    # # Compute SFS vector and CTS vector
    # sfs <- tabulate(colSums(train_mat))
    # cts <- sapply(2:n_train, function(n) ncol(train_mat[1:n,colSums(train_mat[1:n,]) > 0])   )
    # cts <- c(0, sum(train_mat[1,]) , cts)
    # 
    # list_extr_GT[[lab_comb_ibp]] <- predict_good_toulmin(n_train, M, sfs, cts, alternative = 0)$preds

  }
  
  # Save the entire workspace related to the mechanism just performed
  save(list = ls(all.names = TRUE), file =  paste0("R_script_paper/eb_poly_",xi,"_fit_estimate_singledataset.RData"))
}


# Repeated dataset: function to produce fit and estimate on specific ecological scenario ------

fit_estimate_polynomial_scenario_repeateddataset <- function(xi, n_dataset = 10, 
                                                             init_obj_PoissonBB, mcmcparams_obj_PoissonBB, prior_obj_PoissonBB,
                                                             init_obj_NegBinBB, mcmcparams_obj_NegBinBB, prior_obj_NegBinBB, c_fr,
                                                             init_obj_GammaIBP, mcmcparams_obj_GammaIBP, prior_obj_GammaIBP,
                                                             seed = 1234){
  
  # Set training dimensions and dimension of the whole sample 
  Ns <- c(20, 40, 80)
  L <- 500
  # Set maximum number of features
  H <- 10^6
  
  # Structures to save fit and estimates 
  labels_comb_bb <- paste(paste("n_train", Ns, sep = "."),"Nbar.emp", sep=":")
  labels_comb_ibp <- paste("n_train", Ns, sep = ".")
  
  # Dataframe to store the number of hitherto unseen features in the test and the number of distinct in training
  df_obs_new <- data.frame(matrix(nrow = n_dataset, ncol = length(Ns)))
  colnames(df_obs_new) <- labels_comb_ibp
  df_obs_train <- data.frame(matrix(nrow = n_dataset, ncol = length(Ns)))
  colnames(df_obs_train) <- labels_comb_ibp
  
  # Dataframe to store the extrapolation at the last test subject
  df_extr_last_PoissonBB <- data.frame(matrix(nrow = n_dataset, ncol = length(Ns)))
  colnames(df_extr_last_PoissonBB) <- labels_comb_bb
  df_extr_last_NegBinBB <- data.frame(matrix(nrow = n_dataset, ncol = length(Ns)))
  colnames(df_extr_last_NegBinBB) <- labels_comb_bb
  df_extr_last_GammaIBP <- data.frame(matrix(nrow = n_dataset, ncol = length(Ns)))
  colnames(df_extr_last_GammaIBP) <- labels_comb_ibp
  
  # Dataframe to store the Nbar_emp's
  df_Nbar_emp <- data.frame(matrix(nrow = n_dataset, ncol = length(Ns)))
  colnames(df_Nbar_emp) <- labels_comb_bb
  
  # List to store SBSP's extrapolation at the last test subject
  df_extr_last_SBSP <- data.frame(matrix(nrow = n_dataset, ncol = length(Ns)))
  colnames(df_extr_last_SBSP) <- labels_comb_bb
  # List to store smoothed GT's extrapolation at the last test subject
  df_extr_last_GT <- data.frame(matrix(nrow = n_dataset, ncol = length(Ns)))
  colnames(df_extr_last_GT) <- labels_comb_bb
  
  
  # Loop over the different datasets 
  
  for (d in 1:n_dataset){
    
    # Generate data 
    data_mat <- generate_data(xi = xi, n = L, H = H, seed = seed + d)
    
    # Loop over the different training set dimensions 
    
    for (j in 1:length(Ns)){
      
      n_train <- Ns[j]
      M <- L - n_train
      
      train_mat <- data_mat[1:n_train,]
      train_mat <- train_mat[, colSums(is.na(train_mat))==0]
      train_mat <- train_mat[, colSums(train_mat)!=0]
      
      train_list <- convert_features_list(train_mat)
      
      test_mat <- data_mat[(n_train+1):L, ]
      test_mat <- test_mat[, colSums(is.na(test_mat))==0]
      test_mat <- test_mat[, colSums(test_mat)!=0]
      
      test_list <- convert_features_list(test_mat)
      
      # Empirical estimate of E(N) is obtained by Chiu
      Nbar_emp <- beta_binomial_estimator(train_mat)
      df_Nbar_emp[d,j] <- Nbar_emp
      
      # Fit the models
      # PoissonBB
      prior_obj_PoissonBB$lambda <- Nbar_emp
      PoissonBB_fit <- GibbsFA(feature_matrix = train_mat, 
                               model = "PoissonBB", 
                               prior = prior_obj_PoissonBB, 
                               initialization = init_obj_PoissonBB, 
                               mcmcparams = mcmcparams_obj_PoissonBB)
      # NegBinBB
      prior_obj_NegBinBB$n0 <- Nbar_emp/(c_fr - 1)
      prior_obj_NegBinBB$mu0 <- 1/c_fr
      NegBinBB_fit <- GibbsFA(feature_matrix = train_mat, 
                              model = "NegBinBB", 
                              prior = prior_obj_NegBinBB, 
                              initialization = init_obj_NegBinBB, 
                              mcmcparams = mcmcparams_obj_NegBinBB)    
      # GammaIBP
      GammaIBP_fit <- GibbsFA(feature_matrix = train_mat, 
                              model = "GammaIBP", 
                              prior = prior_obj_GammaIBP, 
                              initialization = init_obj_GammaIBP, 
                              mcmcparams = mcmcparams_obj_GammaIBP)  
      
      # Fill the structures about number of features in train and test
      feat_train <- unique(unlist(train_list))
      feat_test <- unique(unlist(test_list))
      df_obs_train[d,j] <- length(feat_train)
      obs_new_features <- setdiff(feat_test, feat_train)
      df_obs_new[d,j] <- length(obs_new_features)
      
      # Fill the extrapolation (last subject) structures
      extr_PoissonBB <- extrapolation(object = PoissonBB_fit, M = M, only_last = TRUE)
      extr_NegBinBB <- extrapolation(object = NegBinBB_fit, M = M, only_last = TRUE)
      extr_GammaIBP <- extrapolation(object = GammaIBP_fit, M = M, only_last = TRUE)
      
      df_extr_last_PoissonBB[d, j] <- mean(unlist(extr_PoissonBB))
      df_extr_last_NegBinBB[d, j] <- mean(unlist(extr_NegBinBB))
      df_extr_last_GammaIBP[d, j] <- mean(unlist(extr_GammaIBP))
      
      
      # Competitors
      
      # A) SBSP's extrapolation
      
      output_sp <- sampler_SB_SP(Z = train_mat,
                                c_0 = 10, beta_0 = 10, alpha_0 = 0.5,
                                p = 0.05, r = 0.1, t = 0.01, a_alpha = 2, b_alpha = 2,
                                tau = 0.005, S = 100, n_burnin = 10, thin = 2, seed = seed)
      
      SBSP_fit <- list("feature_matrix" = train_mat,
                       "alpha_chain" = output_sp$alpha_chain,
                       "theta_chain" = 1 - output_sp$alpha_chain,
                       "a_chain" = output_sp$c_chain +1,
                       "b_chain" = output_sp$beta_chain*(1-output_sp$alpha_chain)/output_sp$alpha_chain)
      
      class(SBSP_fit) <- c("GibbsFA", "GammaIBP")
      
      extr_SBSP <- extrapolation(object = SBSP_fit, M = M, only_last = TRUE)
      
      df_extr_last_SBSP[d, j] <- mean(unlist(extr_SBSP))
      
      
      # B) Smoothed Good-Toulmin extrapolation
      sfs <- tabulate(colSums(train_mat))
      cts <- sapply(2:n_train, function(n) ncol(train_mat[1:n,colSums(train_mat[1:n,]) > 0])   )
      cts <- c(0, sum(train_mat[1,]) , cts)
      
      df_extr_last_GT[d, j] <- predict_good_toulmin(n_train, M, sfs, cts, alternative = 0)$preds [n_train+M+1]
      
      
    }
    
  }
  
  
  # Save the entire workspace related to the mechanism just performed
  save(list = ls(all.names = TRUE), file =  paste0("R_script_paper/poly_",xi,"_fit_estimate_repeateddataset.RData"))
}



########### 1) Main script single-dataset -------

### 1.A) Prior approach ----

# Choose mechanism
xi = 1 # c(0.8, 1.2)

# Fit and estimate richness, rarefaction and extrapolation for GibbsFA's (save workspace)
if (!file.exists(paste0("R_script_paper/poly_",xi,"_fit_estimate_singledataset.RData"))) {
  
  # 1) PoissonBB 
  
  # Initialization and MCMC setting 
  init_PoissonBB <- list(alpha_0 = -1, s_0 = 1)
  init_obj_PoissonBB <- initialization(model = "PoissonBB", init = init_PoissonBB )
  mcmcparams_PoissonBB <- list(tau = 0.1, S = 300, n_burnin = 100, thin = 2)
  mcmcparams_obj_PoissonBB <- mcmcparameters(model = "PoissonBB", mcmcparams = mcmcparams_PoissonBB)
  
  # Hyperparameters elicitation 
  hyper_PoissonBB <- list(a_alpha = 1, b_alpha = 0.1,
                          a_s = 2, b_s = 0.2,
                          lambda = 1)
  prior_obj_PoissonBB <- prior(model = "PoissonBB", hyper = hyper_PoissonBB) 
  
  
  # 2) NegBinBB
  
  # Initialization and MCMC setting
  init_NegBinBB <- list(alpha_0 = -1, s_0 = 1)
  init_obj_NegBinBB <- initialization(model = "NegBinBB", init = init_NegBinBB )
  mcmcparams_NegBinBB <- list(tau = 0.1, S = 300, n_burnin = 100, thin = 2)
  mcmcparams_obj_NegBinBB <- mcmcparameters(model = "NegBinBB", mcmcparams = mcmcparams_NegBinBB)
  
  # Hyperparameters elicitation 
  c_fr <- 10
  
  hyper_NegBinBB <- list(a_alpha = 1, b_alpha = 0.1,
                         a_s = 2, b_s = 0.2,
                         n0 = 1, # n0, mu0 are set s.t. E(N) = Nbar, Var(N) = c_fr*E(N)
                         mu0 = 0.5)
  prior_obj_NegBinBB <- prior(model = "NegBinBB", hyper = hyper_NegBinBB) 
  
  # 3) GammaIBP
  
  # Initialization and MCMC setting 
  init_GammaIBP <- list(alpha_0 = 0.5, s_0 = 15, a_0 = 5, b_0 = 1)
  init_obj_GammaIBP <- initialization(model = "GammaIBP", init = init_GammaIBP )
  mcmcparams_GammaIBP <- list(sigq_alpha = 0.1, sigq_s = 0.1, 
                              S = 300, n_burnin = 100, thin = 2)
  mcmcparams_obj_GammaIBP <- mcmcparameters(model = "GammaIBP", mcmcparams = mcmcparams_GammaIBP)
  
  # Hyperparameters elicitation 
  hyper_GammaIBP <- list(a_alpha = 2, b_alpha = 2,
                         a_s = 2, b_s = 0.2,
                         q = 0.05, r = 1, t = 0.1)
  prior_obj_GammaIBP <- prior(model = "GammaIBP", hyper = hyper_GammaIBP) 
  
  # 4) Call the routine to perform simulations
  
  fit_estimate_polynomial_scenario_singledataset(xi = xi,
                                                 init_obj_PoissonBB, mcmcparams_obj_PoissonBB, prior_obj_PoissonBB,
                                                 init_obj_NegBinBB, mcmcparams_obj_NegBinBB, prior_obj_NegBinBB, c_fr,
                                                 init_obj_GammaIBP, mcmcparams_obj_GammaIBP, prior_obj_GammaIBP )
  
}

# Load the Work space
load(paste0("R_script_paper/poly_",xi,"_fit_estimate_singledataset.RData"))


# 1) Plot Extrapolation

df_extr_PoissonBB_long <- list_extr_GibbsFA_to_long(list_extr_PoissonBB, model = "Poisson")
df_extr_NegBinBB_long <- list_extr_GibbsFA_to_long(list_extr_NegBinBB, model = "NegBin")
df_extr_GammaIBP_long <- list_extr_GibbsFA_to_long(list_extr_GammaIBP, model = "Gamma")

df_extr_GT_long <- list_extr_competitor_to_long(list_extr_GT, model = "GT")  
df_extr_SBSP_long <- list_extr_GibbsFA_to_long(list_extr_SBSP, model = "Gamma") %>%
  mutate(model = "SBSP")

accum <- rarefaction(data_mat)
accum_df <- data.frame("accum" = c(0, accum),
                       "t" = 0:length(accum))

df_GibbsFA <- rbind(df_extr_PoissonBB_long, df_extr_NegBinBB_long, df_extr_GammaIBP_long,
                    df_extr_SBSP_long) %>%
  filter(Nbar %in% c("Not applicable", "emp"))

temp <- tibble(n_train = Ns, xvalues = Ns)

dev.new()
ggplot(df_GibbsFA, aes(x = t, y = means, color = model)) +
  geom_line() +
  geom_line(data = accum_df, aes(t, accum), color="black", linetype="solid") +
  geom_line(data = df_extr_GT_long, aes(t, value)) +
  facet_wrap(.~ n_train, scales = "free_x") +
  geom_vline(data = temp, mapping =  aes(xintercept = xvalues) , linetype = "dashed", color = "grey") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_tableau() +
  theme(aspect.ratio = 1)



### 1.B) EB approach: EFPF vs MM ----

# Choose mechanism
xi = 1 

# Fit and estimate richness, rarefaction and extrapolation for GibbsFA's (save workspace)
if (!file.exists(paste0("R_script_paper/eb_poly_",xi,"_fit_estimate_singledataset.RData"))) {
  
 
  c_fr <- 5
 
  # 3) GammaIBP
  
  # Initialization and known parameters
  eb_init_GammaIBP <- list(alpha = 0.5, s = 1, a = 1, b = 1)
  eb_known_GammaIBP <- list()
  
  # 4) Call the routine to perform simulations
  eb_fit_estimate_polynomial_scenario_singledataset(xi = xi, 
                                                    c_fr,
                                                    eb_init_GammaIBP, eb_known_GammaIBP,
                                                    seed = 123456)
  
}

# Load the Work space
load(paste0("R_script_paper/eb_poly_",xi,"_fit_estimate_singledataset.RData"))

Kn <- sapply(Ns, function(n) sum(colSums(data_mat[1:n,]) > 0)  )
sum(colSums(data_mat) > 0)

# 0) Model checking: verify model is compatible with data 
plot( x = 0:200, y = c(0, rarefaction(data_mat[1:200,])) )

emp_pis <- bind_rows(lapply(Ns, function(n)
  data.frame(x = colMeans(data_mat[1:n,])[colMeans(data_mat[1:n,]) > 0],
             n_train = n)))

# PoissonBB
list_eb_fit_PoissonBB_Nbar_fix <- list_eb_fit_PoissonBB[grepl("Nbar.emp", 
                                                              names(list_eb_fit_PoissonBB))]

params_beta <- data.frame(alpha = sapply(1:length(list_eb_fit_PoissonBB_Nbar_fix), function(i)
  list_eb_fit_PoissonBB_Nbar_fix[[i]]$alpha),
  theta = sapply(1:length(list_eb_fit_PoissonBB_Nbar_fix), function(i)
    list_eb_fit_PoissonBB_Nbar_fix[[i]]$theta),
  lambda = sapply(1:length(list_eb_fit_PoissonBB_Nbar_fix), function(i)
    list_eb_fit_PoissonBB_Nbar_fix[[i]]$lambda ),
  n_train = Ns)

grid <- seq(0,1, length.out = 1000)
# cdf_betas <- bind_rows(lapply(1:nrow(params_beta), function(i)
#   data.frame(x = grid,
#              y = pbeta(grid, shape1 = - params_beta$alpha[i],
#                        shape2 = params_beta$alpha[i] + params_beta$theta[i]),
#              n_train = Ns[i]) ) )

a_beta <- - params_beta$alpha
b_beta <- params_beta$alpha + params_beta$theta

cdf_betas_cond <- bind_rows(lapply(1:nrow(params_beta), function(i)
  data.frame(x = grid,
             y = pbeta(grid, shape1 = a_beta[i], shape2 = b_beta[i])*
               (1/(1- beta(a_beta[i], Ns[i] + b_beta[i])/beta(a_beta[i], b_beta[i]))) +
               pbeta(grid, shape1 = a_beta[i], shape2 = Ns[i] + b_beta[i])*
               (1/(1- beta(a_beta[i], b_beta[i])/beta(a_beta[i], Ns[i] + b_beta[i]))),
             n_train = Ns[i]) ) )


n_train.labs <- paste0("n = ", Ns,", Kn = ", Kn )
names(n_train.labs) <- Ns

ggplot(emp_pis, aes(x = x) ) +
  stat_ecdf(linewidth=2, colour = "red") +
  facet_wrap(.~ n_train, labeller = labeller(n_train = n_train.labs ), scales = "free_x", nrow = 1) +
  geom_line(data = cdf_betas_cond, aes(x = x, y = y), colour = "blue") +
  labs(title="ECDF and theoretical CDF")  


# NegBinBB
list_eb_fit_NegBinBB_Nbar_fix <- list_eb_fit_NegBinBB[grepl("Nbar.emp", 
                                                            names(list_eb_fit_NegBinBB))]

params_beta_NegBinBB <- data.frame(alpha = sapply(1:length(list_eb_fit_NegBinBB_Nbar_fix), function(i)
  list_eb_fit_NegBinBB_Nbar_fix[[i]]$alpha),
  theta = sapply(1:length(list_eb_fit_NegBinBB_Nbar_fix), function(i)
    list_eb_fit_NegBinBB_Nbar_fix[[i]]$theta),
  var_fct = sapply(1:length(list_eb_fit_NegBinBB_Nbar_fix), function(i)
    list_eb_fit_NegBinBB_Nbar_fix[[i]]$var_fct ),
  mu0 = sapply(1:length(list_eb_fit_NegBinBB_Nbar_fix), function(i)
    list_eb_fit_NegBinBB_Nbar_fix[[i]]$mu0 ),
  n_train = Ns)

grid <- seq(0,1, length.out = 1000)
# cdf_betas <- bind_rows(lapply(1:nrow(params_beta), function(i)
#   data.frame(x = grid,
#              y = pbeta(grid, shape1 = - params_beta$alpha[i],
#                        shape2 = params_beta$alpha[i] + params_beta$theta[i]),
#              n_train = Ns[i]) ) )

a_beta <- - params_beta_NegBinBB$alpha
b_beta <- params_beta_NegBinBB$alpha + params_beta_NegBinBB$theta

cdf_betas_cond_NegBinBB <- bind_rows(lapply(1:nrow(params_beta_NegBinBB), function(i)
  data.frame(x = grid,
             y = pbeta(grid, shape1 = a_beta[i], shape2 = b_beta[i])*
               (1/(1- beta(a_beta[i], Ns[i] + b_beta[i])/beta(a_beta[i], b_beta[i]))) +
               pbeta(grid, shape1 = a_beta[i], shape2 = Ns[i] + b_beta[i])*
               (1/(1- beta(a_beta[i], b_beta[i])/beta(a_beta[i], Ns[i] + b_beta[i]))),
             n_train = Ns[i]) ) )


n_train.labs <- paste0("n = ", Ns,", Kn = ", Kn )
names(n_train.labs) <- Ns

ggplot(emp_pis, aes(x = x) ) +
  stat_ecdf(linewidth=2, colour = "red") +
  facet_wrap(.~ n_train, labeller = labeller(n_train = n_train.labs ), scales = "free_x", nrow = 1) +
  geom_line(data = cdf_betas_cond_NegBinBB, aes(x = x, y = y), colour = "blue") +
  labs(title="ECDF and theoretical CDF")  


##### Comparison between PoissonBB(cens), NegBinBB(cens), GammaIBP(EFPF) and GammaIBP(MM) ----

# 0.B) Check on rarefaction
n_rare <- Ns[3]
eb_MM_cens_fit_PoissonBB_rare <- list_eb_MM_cens_fit_PoissonBB[[3]]
eb_MM_cens_fit_NegBinBB_rare <- list_eb_MM_cens_fit_NegBinBB[[3]]
eb_EFPF_fit_GammaIBP_rare <- list_eb_EFPF_fit_GammaIBP[[3]]
eb_MM_fit_GammaIBP_rare <- list_eb_MM_fit_GammaIBP[[3]]

accum_df <- tibble( x = 0:n_rare,
                    n_feat = c(0,rarefaction(data_mat[1:n_rare,], n_reorderings = 10)))

rare_MM_cens_PoissonBB <- tibble( lambda_post = unname(unlist(
  rarefaction(object = eb_MM_cens_fit_PoissonBB_rare, seed = seed)$lambda_post ))) %>%
  rename(means = lambda_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "PoissonBB_MM_cens",
             x = c(1:n_rare,0))

rare_MM_cens_NegBinBB <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_MM_cens_fit_NegBinBB_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "NegBinBB_MM_cens",
             x = c(1:n_rare,0))

rare_EFPF_GammaIBP <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_EFPF_fit_GammaIBP_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "GammaIBP_EFPF",
             x = c(1:n_rare,0))

rare_MM_GammaIBP <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_MM_fit_GammaIBP_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "GammaIBP_MM",
             x = c(1:n_rare,0))

df_rare <- rbind(rare_MM_cens_PoissonBB,
                 rare_MM_cens_NegBinBB,
                 rare_EFPF_GammaIBP, rare_MM_GammaIBP)

df_rare$Model <- factor(df_rare$Model)


ggplot(df_rare, aes(x = x, y = means, color = Model)) +
  geom_line(linetype = "solid", color = "red" , linewidth = 0.9) +
  facet_wrap(.~ Model, scales = "free_x", nrow = 1) +
  #geom_ribbon(aes(ymin = lb_bands, ymax = ub_bands), color = "red" , linewidth = 0.8, alpha = 0.1) +
  geom_point( data = accum_df, aes(x = x, y = n_feat), color="black", shape = 1, size = 0.5) +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1)
ggsave(filename = paste0("R_script_paper/Paper_plots/rarefaction_poly_", xi, "_eb.pdf"), width = 8, height = 3.5, dpi = 300, units = "in", device='pdf')


# 0.c) Check on K_n_r
n_knr <- Ns[3]
eb_MM_cens_fit_PoissonBB_knr <- list_eb_MM_cens_fit_PoissonBB[[3]]
eb_MM_cens_fit_NegBinBB_knr <- list_eb_MM_cens_fit_NegBinBB[[3]]
eb_EFPF_fit_GammaIBP_knr <- list_eb_EFPF_fit_GammaIBP[[3]]
eb_MM_fit_GammaIBP_knr <- list_eb_MM_fit_GammaIBP[[3]]


observed_K_n_r <- tibble( r = 1:n_knr,
                          k_n_r = K_n_r(data_mat[1:n_knr,], n_reorderings = 1)[[paste0('N = ', n_knr)]])

K_n_r_MM_cens_PoissonBB <- tibble( lambda_est = unname(unlist(
  K_n_r(object = eb_MM_cens_fit_PoissonBB_knr, n = n_knr)[[paste0('N = ', n_knr)]]$lambda_est ))) %>%
  rename(means = lambda_est) %>%
  add_column(Model = "PoissonBB_MM_cens",
             r = 1:n_knr)

K_n_r_MM_cens_NegBinBB <- tibble( mu0_est = unname(unlist(
  K_n_r(object = eb_MM_cens_fit_NegBinBB_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
  rename(means = mu0_est) %>%
  add_column(Model = "NegBinBB_MM_cens",
             r = 1:n_knr)

K_n_r_EFPF_GammaIBP <- tibble( mu0_est = unname(unlist(
  K_n_r(object = eb_EFPF_fit_GammaIBP_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
  rename(means = mu0_est) %>%
  add_column(Model = "GammaIBP_EFPF",
             r = 1:n_knr)

K_n_r_MM_GammaIBP <- tibble( mu0_est = unname(unlist(
  K_n_r(object = eb_MM_fit_GammaIBP_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
  rename(means = mu0_est) %>%
  add_column(Model = "GammaIBP_MM",
             r = 1:n_knr)

df_K_n_r <- rbind(K_n_r_MM_cens_PoissonBB,
                  #K_n_r_MM_cens_NegBinBB,
                  K_n_r_EFPF_GammaIBP, K_n_r_MM_GammaIBP)

df_K_n_r$Model <- factor(df_K_n_r$Model)

r_positive <- observed_K_n_r %>%
  filter(k_n_r > 0) %>%
  select(r) %>%
  filter(r < 10)

df_K_n_r_plot <- df_K_n_r %>%
  filter(r %in% c(r_positive$r))

observed_K_n_r_plot <- observed_K_n_r %>%
  filter(r %in% c(r_positive$r))


ggplot(df_K_n_r_plot, aes(x = r, y = means, color = Model)) +
  geom_line( linetype = "dashed") +
  geom_point( data = observed_K_n_r_plot, aes(x = r, y = k_n_r), color="black", size = 1) +
  scale_y_log10() +
  xlab("r") + ylab(expression(m[r])) + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1)
ggsave(filename = paste0("R_script_paper/Paper_plots/knr_poly_", xi, "_eb.pdf"), width = 5, height = 4.5, dpi = 300, units = "in", device='pdf')




# 3) Plot Extrapolation - EB version

# Extract accumulation curve of the observed sample (or average accumulation)
M = 400

accum_df <- tibble(n_feat = unlist(sapply(Ns, function(n) 
  c(0,rarefaction(data_mat[1:(n + M), ], n_reorderings = 1)))),
  n_train = rep(Ns, times = Ns + M +1),
  type = unlist(sapply(Ns, function(n) c(rep("train", n +1), rep("test", M)))),
  x = unlist(sapply(Ns, function(n) 0:(n + M))))

accum_df_train <- accum_df %>%
  filter(type == "train")
accum_df_test <- accum_df %>%
  filter(type == "test")


# Poisson
extr_MM_cens_PoissonBB_df <- tibble(lambda = unname(unlist(lapply(list_eb_MM_cens_PoissonBB, function(x)
  extrapolation(object = x, M = M, seed = seed)$lambda_post ))),
  n_train = rep(Ns, each = M),
  Kn = rep(Kn, each = M)) %>%
  mutate(lb = qpois(0.025, lambda, lower.tail = TRUE, log.p = FALSE),
         ub = qpois(0.975, lambda, lower.tail = TRUE, log.p = FALSE)) %>%
  rename(means = lambda) %>%
  add_row(means = 0, lb = 0, ub = 0, n_train = Ns, Kn = Kn) %>%
  mutate(means = means + Kn, lb = lb + Kn, ub = ub + Kn) %>%
  add_column(x = c(unlist(sapply(Ns, function(n) (n+1):(M+n))), Ns),
             Model = "PoissonBB_MM_cens")

extr_MM_cens_PoissonBB_df$x <- as.integer(extr_MM_cens_PoissonBB_df$x)
extr_MM_cens_PoissonBB_df <- extr_MM_cens_PoissonBB_df %>%
  select(means, lb, ub, n_train, x, Model)

# NegBin
extr_MM_cens_NegBinBB_df <- tibble(mu0 = unname(unlist(lapply(list_eb_MM_cens_NegBinBB, function(x)
  extrapolation(object = x, M = M, seed = seed)$mu0_post ))),
  n0 = unname(unlist(lapply(list_eb_MM_cens_NegBinBB, function(x)
    extrapolation(object = x, M = M, seed = seed)$n0_post ))),
  n_train = rep(Ns, each = M),
  Kn = rep(Kn, each = M)) %>%
  mutate(p = 1/(mu0/n0 + 1),
         lb = qnbinom(0.025, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE),
         ub = qnbinom(0.975, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE)) %>%
  rename(means = mu0) %>%
  add_row(means = 0, lb = 0, ub = 0, n_train = Ns, Kn = Kn) %>%
  mutate(means = means + Kn, lb = lb + Kn, ub = ub + Kn) %>%
  add_column(x = c(unlist(sapply(Ns, function(n) (n+1):(M+n))), Ns),
             Model = "NegBinBB_MM_cens")

extr_MM_cens_NegBinBB_df$x <- as.integer(extr_MM_cens_NegBinBB_df$x)
extr_MM_cens_NegBinBB_df <- extr_MM_cens_NegBinBB_df %>%
  select(means, lb, ub, n_train, x, Model)


# GammaIBP
extr_EFPF_GammaIBP_df <- tibble(mu0 = unname(unlist(lapply(list_eb_EFPF_GammaIBP, function(x)
  extrapolation(object = x, M = M, seed = seed)$mu0_post ))),
  n0 = unname(unlist(lapply(list_eb_EFPF_GammaIBP, function(x)
    extrapolation(object = x, M = M, seed = seed)$n0_post ))),
  n_train = rep(Ns, each = M),
  Kn = rep(Kn, each = M)) %>%
  mutate(p = 1/(mu0/n0 + 1),
         lb = qnbinom(0.025, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE),
         ub = qnbinom(0.975, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE)) %>%
  rename(means = mu0) %>%
  add_row(means = 0, lb = 0, ub = 0, n_train = Ns, Kn = Kn) %>%
  mutate(means = means + Kn, lb = lb + Kn, ub = ub + Kn) %>%
  add_column(x = c(unlist(sapply(Ns, function(n) (n+1):(M+n))), Ns),
             Model = "GammaIBP_EFPF")

extr_EFPF_GammaIBP_df$x <- as.integer(extr_EFPF_GammaIBP_df$x)
extr_EFPF_GammaIBP_df <- extr_EFPF_GammaIBP_df %>%
  select(means, lb, ub, n_train, x, Model)


extr_MM_GammaIBP_df <- tibble(mu0 = unname(unlist(lapply(list_eb_MM_GammaIBP, function(x)
  extrapolation(object = x, M = M, seed = seed)$mu0_post ))),
  n0 = unname(unlist(lapply(list_eb_MM_GammaIBP, function(x)
    extrapolation(object = x, M = M, seed = seed)$n0_post ))),
  n_train = rep(Ns, each = M),
  Kn = rep(Kn, each = M)) %>%
  mutate(p = 1/(mu0/n0 + 1),
         lb = qnbinom(0.025, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE),
         ub = qnbinom(0.975, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE)) %>%
  rename(means = mu0) %>%
  add_row(means = 0, lb = 0, ub = 0, n_train = Ns, Kn = Kn) %>%
  mutate(means = means + Kn, lb = lb + Kn, ub = ub + Kn) %>%
  add_column(x = c(unlist(sapply(Ns, function(n) (n+1):(M+n))), Ns),
             Model = "GammaIBP_MM")

extr_MM_GammaIBP_df$x <- as.integer(extr_MM_GammaIBP_df$x)
extr_MM_GammaIBP_df <- extr_MM_GammaIBP_df %>%
  select(means, lb, ub, n_train, x, Model)


# Plot
n_train.labs <- paste0("n = ", Ns,", Kn = ", Kn )
names(n_train.labs) <- Ns
temp <- tibble(n_train = Ns, xvalues = Ns)
extr_all_df <- rbind(extr_PoissonBB_df, extr_NegBinBB_df, extr_GammaIBP_df)

ggplot(extr_all_df, aes(x, means, color = Model)) +
  geom_line(linetype = "dashed", linewidth = 0.9) +
  geom_point( data = accum_df_train, aes(x = x, y = n_feat),
              color="black", shape = 1, size = 1) +
  geom_point( data = accum_df_test, aes(x = x, y = n_feat),
              color="darkgrey", shape = 1, size = 0.1) +
  #geom_line(data = df_extr_GT_long, aes(t, value), linetype = "dashed", linewidth = 1) +
  facet_wrap(. ~ n_train,labeller = labeller(n_train = n_train.labs ),  scales = "free_x") +
  geom_vline(data = temp, mapping =  aes(xintercept = xvalues) , linetype = "dashed", color = "grey") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1)
ggsave(filename = "R_script_paper/Paper_plots/extr_poly_1_eb.pdf", width = 8, height = 3.8, dpi = 300, units = "in", device='pdf')



# 2) Main script repeated-dataset -------

# Choose mechanism
xi = 1 # c(0.8, 1.2)
n_dataset = 2

# Fit and estimate richness, rarefaction and extrapolation for GibbsFA's (save workspace)
if (!file.exists(paste0("R_script_paper/poly_",xi,"_fit_estimate_repeateddataset.RData"))) {
  
  # 1) PoissonBB 
  
  # Initialization and MCMC setting 
  init_PoissonBB <- list(alpha_0 = -1, s_0 = 1)
  init_obj_PoissonBB <- initialization(model = "PoissonBB", init = init_PoissonBB )
  mcmcparams_PoissonBB <- list(tau = 0.1, S = 300, n_burnin = 100, thin = 2)
  mcmcparams_obj_PoissonBB <- mcmcparameters(model = "PoissonBB", mcmcparams = mcmcparams_PoissonBB)
  
  # Hyperparameters elicitation 
  hyper_PoissonBB <- list(a_alpha = 1, b_alpha = 0.1,
                          a_s = 2, b_s = 0.2,
                          lambda = 1)
  prior_obj_PoissonBB <- prior(model = "PoissonBB", hyper = hyper_PoissonBB) 
  
  
  # 2) NegBinBB
  
  # Initialization and MCMC setting
  init_NegBinBB <- list(alpha_0 = -1, s_0 = 1)
  init_obj_NegBinBB <- initialization(model = "NegBinBB", init = init_NegBinBB )
  mcmcparams_NegBinBB <- list(tau = 0.1, S = 300, n_burnin = 100, thin = 2)
  mcmcparams_obj_NegBinBB <- mcmcparameters(model = "NegBinBB", mcmcparams = mcmcparams_NegBinBB)
  
  # Hyperparameters elicitation 
  c_fr <- 10
  
  hyper_NegBinBB <- list(a_alpha = 1, b_alpha = 0.1,
                         a_s = 2, b_s = 0.2,
                         n0 = 1, # n0, mu0 are set s.t. E(N) = Nbar, Var(N) = c_fr*E(N)
                         mu0 = 0.5)
  prior_obj_NegBinBB <- prior(model = "NegBinBB", hyper = hyper_NegBinBB) 
  
  # 3) GammaIBP
  
  # Initialization and MCMC setting 
  init_GammaIBP <- list(alpha_0 = 0.5, s_0 = 15, a_0 = 5, b_0 = 1)
  init_obj_GammaIBP <- initialization(model = "GammaIBP", init = init_GammaIBP )
  mcmcparams_GammaIBP <- list(sigq_alpha = 0.1, sigq_s = 0.1, 
                              S = 300, n_burnin = 100, thin = 2)
  mcmcparams_obj_GammaIBP <- mcmcparameters(model = "GammaIBP", mcmcparams = mcmcparams_GammaIBP)
  
  # Hyperparameters elicitation 
  hyper_GammaIBP <- list(a_alpha = 2, b_alpha = 2,
                         a_s = 2, b_s = 0.2,
                         q = 0.05, r = 1, t = 0.1)
  prior_obj_GammaIBP <- prior(model = "GammaIBP", hyper = hyper_GammaIBP) 
  
  # 4) Call the routine to perform simulations
  fit_estimate_polynomial_scenario_repeateddataset(xi = xi, n_dataset = n_dataset,
                                                   init_obj_PoissonBB, mcmcparams_obj_PoissonBB, prior_obj_PoissonBB,
                                                   init_obj_NegBinBB, mcmcparams_obj_NegBinBB, prior_obj_NegBinBB, c_fr,
                                                   init_obj_GammaIBP, mcmcparams_obj_GammaIBP, prior_obj_GammaIBP)
}

# Load the Work space
load(paste0("R_script_paper/poly_",xi,"_fit_estimate_repeateddataset.RData"))



# 1) Plot boxplots on accuracy 

df_new_PoissonBB <- df_extr_last_PoissonBB - df_obs_train
df_new_NegBinBB <- df_extr_last_NegBinBB - df_obs_train
df_new_GammaIBP <- df_extr_last_GammaIBP - df_obs_train
df_new_GT <- df_extr_last_GT - df_obs_train
df_new_SBSP <- df_extr_last_SBSP - df_obs_train

# compute accuracy
acc_alt_PoissonBB <- compute_accuracy(df_obs_new, df_new_PoissonBB, df_obs_train)
acc_alt_NegBinBB <- compute_accuracy(df_obs_new, df_new_NegBinBB, df_obs_train)
acc_alt_GammaIBP <- compute_accuracy(df_obs_new, df_new_GammaIBP, df_obs_train)
acc_alt_GT <- compute_accuracy(df_obs_new, df_new_GT, df_obs_train)
acc_alt_SBSP <- compute_accuracy(df_obs_new, df_new_SBSP, df_obs_train)

# handle the dataframes
acc_alt_PoissonBB_long <- acc_alt_PoissonBB %>%
  pivot_longer( everything(), names_to = "training", values_to = "Accuracy") %>%
  add_column(model = "Poisson", n_train = rep(Ns, each = n_dataset))

acc_alt_NegBinBB_long <- acc_alt_NegBinBB %>%
  pivot_longer( everything(), names_to = "training", values_to = "Accuracy") %>%
  add_column(model = "NegBin", n_train = rep(Ns, each = n_dataset))

acc_alt_GammaIBP_long <- acc_alt_GammaIBP %>%
  pivot_longer( everything(), names_to = "training", values_to = "Accuracy") %>%
  add_column(model = "Gamma", n_train = rep(Ns, each = n_dataset))

acc_alt_GT_long <- acc_alt_GT %>%
  pivot_longer( everything(), names_to = "training", values_to = "Accuracy") %>%
  add_column(model = "GT", n_train = rep(Ns, each = n_dataset))

acc_alt_SBSP_long <- acc_alt_SBSP %>%
  pivot_longer( everything(), names_to = "training", values_to = "Accuracy") %>%
  add_column(model = "SBSP", n_train = rep(Ns, each = n_dataset))


joint_alt_long <- bind_rows(acc_alt_PoissonBB_long, acc_alt_NegBinBB_long,
                            acc_alt_GammaIBP_long, acc_alt_GT_long, acc_alt_SBSP_long)

# plots
ggplot(joint_alt_long, aes(x = model, y=Accuracy)) +
  geom_boxplot() +
  facet_wrap(.~n_train, scales = "free_x", nrow = 1) +
  theme_light() +
  rremove("xlab") +
  ylab("Error index") +
  scale_y_continuous(breaks = pretty_breaks())+
  theme(aspect.ratio = 1)


