
rm(list=ls())
library(ProductFormFA)
library(tidyverse)
library(scales)
library(ggpubr)
library(ggthemes)

source("R_script_paper/Routine_Chao.R")
source("R_script_paper/utils.R")


# Function to generate binary matrix according to 4 different mechanisms ----

generate_data <- function(mechanism, n, H = 100, seed = 1234){
  
  if (!mechanism %in% c("custom", "beta_pis", "homogeneous", "random_uniform", "broken_stick", "log-normal")){
    stop("Invalid generating mechanism.")
  }
  
  set.seed(seed)
  
  if (mechanism == "custom"){
    mul <- 1
    data_mat <- cbind(
      
      # matrix(rbinom(n * 200*mul, size = 1, prob = 0.02), n, 200*mul),
      # matrix(rbinom(n * 200*mul, size = 1, prob = 0.01), n, 200*mul),
      # matrix(rbinom(n * 200*mul, size = 1, prob = 0.005), n, 200*mul)
      matrix(rbinom(n * 200*mul, size = 1, prob = 0.015), n, 200*mul),
      matrix(rbinom(n * 200*mul, size = 1, prob = 0.01), n, 200*mul),
      matrix(rbinom(n * 200*mul, size = 1, prob = 0.005), n, 200*mul)
      
    )
  
    # This is good for GammaIBP
    
    # data_mat <- cbind(
    #   matrix(rbinom(n * 5, size = 1, prob = 0.8), n, 5),
    #   #matrix(rbinom(n * 5, size = 1, prob = 0.7), n, 5),
    #   matrix(rbinom(n * 10, size = 1, prob = 0.6), n, 10),
    #   #matrix(rbinom(n * 10, size = 1, prob = 0.5), n, 10),
    #   matrix(rbinom(n * 10, size = 1, prob = 0.4), n, 10),
    #   #matrix(rbinom(n * 10, size = 1, prob = 0.3), n, 10),
    #   matrix(rbinom(n * 10, size = 1, prob = 0.2), n, 10),
    #   #matrix(rbinom(n * 10, size = 1, prob = 0.1), n, 10),
    #   matrix(rbinom(n * 40, size = 1, prob = 0.05), n, 40),
    #   matrix(rbinom(n * 80, size = 1, prob = 0.01), n, 80),
    #   matrix(rbinom(n * 200, size = 1, prob = 0.001), n, 200)
    # )
    

  } else if (mechanism == "beta_pis"){
    
    pis <- rbeta(H, 1, 100) # 1,20 is good
    
    data_mat <- matrix(rbinom(n*H, size = 1, prob = rep(pis, n)),
                       nrow = n, ncol = H, byrow = T )
    
  } else if (mechanism == "homogeneous"){ 
    
    pres_prob <- 0.05
    
    data_mat <- matrix(rbinom(n*H, size = 1, prob = pres_prob), nrow = n, ncol = H)
    
  } else if (mechanism == "random_uniform"){
    
    as <- runif(H, 0, 1)
    c <- 0.5 / max(as)
    pis <- c*as
    
    data_mat <- matrix(rbinom(n*H, size = 1, prob = rep(pis, n)),
                       nrow = n, ncol = H, byrow = T )
    
  } else if (mechanism == "broken_stick"){
    
    as <- rexp(H, 1)
    c <- 0.5 / max(as)
    pis <- c*as
    
    data_mat <- matrix(rbinom(n*H, size = 1, prob = rep(pis, n)),
                       nrow = n, ncol = H, byrow = T )
    
  } else if (mechanism == "log-normal"){
    
    as <- exp(rnorm(H, 0, 1))
    c <- 1 / max(as)
    pis <- c*as
    
    data_mat <- matrix(rbinom(n*H, size = 1, prob = rep(pis, n)),
                       nrow = n, ncol = H, byrow = T )
    
  }
  
  return(data_mat)
  
}



# Single dataset: function to produce fit and estimate on specific ecological scenario ------

fit_estimate_ecological_scenario_singledataset <- function(mechanism, 
                                                           init_obj_PoissonBB, mcmcparams_obj_PoissonBB, prior_obj_PoissonBB,
                                                           init_obj_NegBinBB, mcmcparams_obj_NegBinBB, prior_obj_NegBinBB, c_fr,
                                                           init_obj_GammaIBP, mcmcparams_obj_GammaIBP, prior_obj_GammaIBP,
                                                           seed = 1234){
  
  if (!mechanism %in% c("homogeneous", "random_uniform", "broken_stick", "log-normal")){
    stop("Invalid generating mechanism.")
  }
  
  # Set training dimensions and dimension of the whole sample 
  Ns <- c(20, 40, 80)
  if (mechanism == "log-normal"){
    L <- 300
  } else {
    L <- 150
  }
  # Set maximum number of features
  H <- 500
  
  # Generate data 
  data_mat <- generate_data(mechanism = mechanism, n = L, H = H, seed = seed)
  
  # Set grid of desired value of E[N] = Nbar, with the empirical case to be added
  Nbars <- c(200, 400, 600)
  
  # Structures to save fit and estimates 
  
  # List to store the MCMC chains for PoissonBB, NegBinBB and GammaIBP, under different settings 
  list_params_PoissonBB <- vector(mode = "list")
  list_params_NegBinBB <-  vector(mode = "list")
  list_params_GammaIBP <-  vector(mode = "list")
  
  # List to store richness draws for PoissonBB and NegBinBB, for different settings
  list_richness_PoissonBB <- vector(mode = "list")
  list_richness_NegBinBB <-  vector(mode = "list")
  
  # List to store rarefaction draws for PoissonBB, NegBinBB and GammaIBP, for different settings
  list_rare_PoissonBB <- vector(mode = "list")
  list_rare_NegBinBB <-  vector(mode = "list")
  list_rare_GammaIBP <-  vector(mode = "list")
  
  # List to store extrapolation draws for PoissonBB, NegBinBB and GammaIBP, for different settings
  list_extr_PoissonBB <- vector(mode = "list")
  list_extr_NegBinBB <-  vector(mode = "list")
  list_extr_GammaIBP <-  vector(mode = "list")
  
  # List of Nbar_emp for the different training sets
  list_Nbar_emp <- vector(mode = "list")
  
  # List to store Chao's estimates
  list_rare_extr_Chao <- vector(mode="list")
  # List to store smoothed GT's estimates
  list_extr_GT <- vector(mode="list")
  
  
  
  # Loop over the different training set dimensions 
  
  for (j in 1:length(Ns)){
    
    n_train <- Ns[j]
    M <- L - n_train
    
    train_mat <- data_mat[1:n_train,]
    
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
    prior_obj_NegBinBB$mu0 <- Nbar_emp
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
    
    # Fill the richness structures
    list_richness_PoissonBB[[lab_comb_bb]] <- total_richness(object = PoissonBB_fit)
    list_richness_NegBinBB[[lab_comb_bb]] <- total_richness(object = NegBinBB_fit)
    
    # Fill the rarefaction structures
    list_rare_PoissonBB[[lab_comb_bb]] <- rarefaction(object = PoissonBB_fit, seed = seed)
    list_rare_NegBinBB[[lab_comb_bb]] <- rarefaction(object = NegBinBB_fit, seed = seed)
    list_rare_GammaIBP[[lab_comb_ibp]] <- rarefaction(object = GammaIBP_fit, seed = seed)
    
    # Fill the extrapolation structures
    list_extr_PoissonBB[[lab_comb_bb]] <- extrapolation(object = PoissonBB_fit, M = M, seed = seed)
    list_extr_NegBinBB[[lab_comb_bb]] <- extrapolation(object = NegBinBB_fit, M = M, seed = seed)
    list_extr_GammaIBP[[lab_comb_ibp]] <- extrapolation(object = GammaIBP_fit, M = M, seed = seed)
    
    
    # Loop over the different Nbar hyperparameters (only for BB mixtures) 
    
    for (v in 1:length(Nbars)){
      
      Nbar <- Nbars[v]
      
      # Fit the models
      # PoissonBB
      prior_obj_PoissonBB$lambda <- Nbar
      PoissonBB_fit <- GibbsFA(feature_matrix = train_mat, 
                               model = "PoissonBB", 
                               prior = prior_obj_PoissonBB, 
                               initialization = init_obj_PoissonBB, 
                               mcmcparams = mcmcparams_obj_PoissonBB)
      # NegBinBB
      prior_obj_NegBinBB$n0 <- Nbar/(c_fr - 1)
      prior_obj_NegBinBB$mu0 <- Nbar
      NegBinBB_fit <- GibbsFA(feature_matrix = train_mat, 
                              model = "NegBinBB", 
                              prior = prior_obj_NegBinBB, 
                              initialization = init_obj_NegBinBB, 
                              mcmcparams = mcmcparams_obj_NegBinBB)  
      
      # Fill the structures of MCMC chains of the parameters
      lab_comb_bb <- paste0("n_train.",n_train,":Nbar.", Nbar)
      
      list_params_PoissonBB[[lab_comb_bb]] <- list("alpha" = PoissonBB_fit$alpha_chain, 
                                                   "theta" = PoissonBB_fit$theta_chain)
      
      list_params_NegBinBB[[lab_comb_bb]] <- list("alpha" = NegBinBB_fit$alpha_chain, 
                                                  "theta" = NegBinBB_fit$theta_chain)
      
      # Fill the richness structures
      list_richness_PoissonBB[[lab_comb_bb]] <- total_richness(object = PoissonBB_fit)
      list_richness_NegBinBB[[lab_comb_bb]] <- total_richness(object = NegBinBB_fit)
      
      # Fill the rarefaction structures
      list_rare_PoissonBB[[lab_comb_bb]] <- rarefaction(object = PoissonBB_fit)
      list_rare_NegBinBB[[lab_comb_bb]] <- rarefaction(object = NegBinBB_fit)
      
      # Fill the extrapolation structures
      list_extr_PoissonBB[[lab_comb_bb]] <- extrapolation(object = PoissonBB_fit, M = M)
      list_extr_NegBinBB[[lab_comb_bb]] <- extrapolation(object = NegBinBB_fit, M = M)
      
      
    }
    
    # Competitors 
    
    # A) Chao's rarefaction and extrapolation
    
    # Determine the frequency vector of the training sets
    Q_vec <- colSums(train_mat)
    Q_vec <- Q_vec[Q_vec>0]
    
    # Compute the curves with confidence intervals
    fit_Chao <- iNEXT.Sam(Spec = Q_vec, T = n_train, endpoint = L)
    
    rare_Chao <- as_tibble(fit_Chao[["q=0"]]) %>%
      select(-Cov.hat) %>%
      rename(medians = D0.hat, lbs = Norm.CI.Low, ubs = Norm.CI.High)
    
    list_rare_extr_Chao[[lab_comb_ibp]] <- as.data.frame(rare_Chao)
    
    
    # B) Smoothed Good-Toulmin extrapolation
    
    # Compute SFS vector and CTS vector
    sfs <- tabulate(colSums(train_mat))
    cts <- sapply(2:n_train, function(n) ncol(train_mat[1:n,colSums(train_mat[1:n,]) > 0])   )
    cts <- c(0, sum(train_mat[1,]) , cts)
    
    list_extr_GT[[lab_comb_ibp]] <- predict_good_toulmin(n_train, M, sfs, cts, alternative = 0)$preds
    
  }
  
  # Save the entire workspace related to the mechanism just performed
  save(list = ls(all.names = TRUE), file =  paste0("R_script_paper/",mechanism,"_fit_estimate_singledataset.RData"))
}




eb_EFPF_fit_estimate_ecological_scenario_singledataset <- function(mechanism,
                                                                   vars_fct_NegBinBB,
                                                                   vars_GammaIBP,
                                                                   eb_init_BB, eb_known_BB,
                                                                   eb_init_IBP, eb_known_IBP,
                                                                   seed = 1234){
  
  if (!mechanism %in% c("custom", "beta_pis", "homogeneous", "random_uniform", "broken_stick", "log-normal")){
    stop("Invalid generating mechanism.")
  }
  
  # Set training dimensions and dimension of the whole sample 
  if (mechanism == "custom"){
    Ns <- c(50, 250, 1250)
  } else {
    Ns <- c(200, 1000, 5000) 
  }
  
  N_max <- Ns[length(Ns)]
  
  if (mechanism == "log-normal"){
    L <- 600
  } else {
    L <- N_max + 600
  }
  # Set maximum number of features
  H <- 500
  
  # Generate data 
  data_mat <- generate_data(mechanism = mechanism, n = L, H = H, seed = seed)
  if (mechanism == "custom"){
    H <- ncol(data_mat) # total features, counting absent ones
  }
  
  data_mat <- data_mat[, colSums(data_mat) > 0]
    
  # Set grid of desired value of E[N] = Nbar, with the empirical case to be added
  if (mechanism == "custom"){
    H_hundred <- round(H / 100) * 100
    Nbars <- c(H_hundred - 400, H_hundred - 200, H_hundred + 200)
  }
  
  # Structures to save fit and estimates 
  
  # List to store the fitted models for PoissonBB, NegBinBB and GammaIBP, under different settings 
  list_eb_EFPF_fit_PoissonBB <- vector(mode = "list")
  
  list_eb_EFPF_fit_NegBinBB <-  vector(mode = "list", length = length(vars_fct_NegBinBB))
  names(list_eb_EFPF_fit_NegBinBB) <- paste0("var_fct.", vars_fct_NegBinBB)
  
  list_eb_EFPF_fit_GammaIBP <-  vector(mode = "list", length = length(vars_GammaIBP))
  names(list_eb_EFPF_fit_GammaIBP) <- paste0("var.", vars_GammaIBP)
  
  
  # Loop over the different training set dimensions 
  
  for (j in 1:length(Ns)){
    
    n_train <- Ns[j]
    M <- L - n_train
    
    train_mat <- data_mat[1:n_train,]
    train_mat <- train_mat[, colSums(train_mat) > 0]
    K <- ncol(train_mat)
    
    lab_comb_bb <- paste0("n_train.",n_train,":Nbar.emp")
    lab_comb_ibp <- paste("n_train", n_train, sep = ".")
    
    eb_params_obj_BB <- eb_params(model = "BB", 
                                  init = eb_init_BB, known = eb_known_BB )
    eb_params_obj_IBP <- eb_params(model = "IBP", 
                                  init = eb_init_IBP, known = eb_known_IBP )
    
    # Fit the models
    # PoissonBB
    list_eb_EFPF_fit_PoissonBB[[lab_comb_bb]] <- GibbsFA_eb(feature_matrix = train_mat, 
                                                            model = "PoissonBB", 
                                                            type = "EFPF",
                                                            eb_params =  eb_params_obj_BB)
    
    # NegBinBB
    for (var_fct_NegBinBB in vars_fct_NegBinBB){
      
      list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]][[lab_comb_bb]] <- 
        GibbsFA_eb(feature_matrix = train_mat,
                   model = "NegBinBB", type = "EFPF",
                   eb_params =  eb_params_obj_BB, 
                   var_fct = var_fct_NegBinBB)
      
    }
    
    
    # GammaIBP
    for (var_GammaIBP in vars_GammaIBP){
      
      list_eb_EFPF_fit_GammaIBP[[paste0("var.", var_GammaIBP)]][[lab_comb_ibp]] <-
        GibbsFA_eb(feature_matrix = train_mat,
                   model = "GammaIBP", type = "EFPF",
                   eb_params =  eb_params_obj_IBP,
                   var_GammaIBP = var_GammaIBP)
      
    }
    
    
    # Loop over the different Nbar hyperparameters (only for BB mixtures) 
    
    for (v in 1:length(Nbars)){
      
      Nbar <- Nbars[v]
      lab_comb_bb <- paste0("n_train.",n_train,":Nbar.", Nbar)
      
      # Fit the models
      # PoissonBB
      eb_known_PoissonBB_Nbar <- list("lambda" = Nbar)
      eb_init_PoissonBB_Nbar <- eb_init_BB[! names(eb_init_BB) %in% c("Nhat_prime")]
      eb_params_obj_PoissonBB <- eb_params(model = "PoissonBB",
                                           init = eb_init_PoissonBB_Nbar,
                                           known = eb_known_PoissonBB_Nbar)
      
      
      list_eb_EFPF_fit_PoissonBB[[lab_comb_bb]] <- 
        GibbsFA_eb(feature_matrix = train_mat, 
                   model = "PoissonBB", 
                   type = "EFPF",
                   eb_params =  eb_params_obj_PoissonBB)
      
      # NegBinBB
      for (var_fct_NegBinBB in vars_fct_NegBinBB){
        
        eb_known_NegBinBB_Nbar <- list("mu0" = Nbar, "var_fct" = var_fct_NegBinBB)
        eb_init_NegBinBB_Nbar <- eb_init_BB[! names(eb_init_BB) %in% c("Nhat_prime")]
        eb_params_obj_NegBinBB <- eb_params(model = "NegBinBB",
                                            init = eb_init_NegBinBB_Nbar,
                                            known = eb_known_NegBinBB_Nbar)
        
        list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]][[lab_comb_bb]] <-
          GibbsFA_eb(feature_matrix = train_mat,
                     model = "NegBinBB", type = "EFPF",
                     eb_params =  eb_params_obj_NegBinBB,
                     var_fct = var_fct_NegBinBB)
        
      }
        
      
                                    
    }
    
    
  }
  
  # Save the entire workspace related to the mechanism just performed
  save(list = ls(all.names = TRUE), file =  paste0("R_script_paper/eb_EFPF_",mechanism,"_fit_estimate_singledataset.RData"))
}


# Repeated dataset: function to produce fit and estimate on specific ecological scenario ------

fit_estimate_ecological_scenario_repeateddataset <- function(mechanism, n_dataset = 10, 
                                                             init_obj_PoissonBB, mcmcparams_obj_PoissonBB, prior_obj_PoissonBB,
                                                             init_obj_NegBinBB, mcmcparams_obj_NegBinBB, prior_obj_NegBinBB, c_fr,
                                                             init_obj_GammaIBP, mcmcparams_obj_GammaIBP, prior_obj_GammaIBP,
                                                             seed = 1234){
  
  if (!mechanism %in% c("homogeneous", "random_uniform", "broken_stick", "log-normal")){
    stop("Invalid generating mechanism.")
  }
  
  # Set training dimensions and dimension of the whole sample 
  Ns <- c(20, 40, 80)
  L <- 150
  # Set maximum number of features
  H <- 500
  
  # Structures to save fit and estimates 
  labels_comb_bb <- paste(paste("n_train", Ns, sep = "."),"Nbar.emp", sep=":")
  labels_comb_ibp <- paste("n_train", Ns, sep = ".")
  
  # Dataframe to store the E(richness) for PoissonBB and NegBinBB, under different settings ad datasets
  df_ev_richness_PoissonBB <- data.frame(matrix(nrow = n_dataset, ncol = length(Ns)))
  colnames(df_ev_richness_PoissonBB) <- labels_comb_bb
  df_ev_richness_NegBinBB <- data.frame(matrix(nrow = n_dataset, ncol = length(Ns)))
  colnames(df_ev_richness_NegBinBB) <- labels_comb_bb
  
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
  
  # List to store Chao's extrapolation at the last test subject
  df_extr_last_Chao <- data.frame(matrix(nrow = n_dataset, ncol = length(Ns)))
  colnames(df_extr_last_Chao) <- labels_comb_bb
  # List to store smoothed GT's extrapolation at the last test subject
  df_extr_last_GT <- data.frame(matrix(nrow = n_dataset, ncol = length(Ns)))
  colnames(df_extr_last_GT) <- labels_comb_bb
  
  
  # Loop over the different datasets 
  
  for (d in 1:n_dataset){
    
    # Generate data 
    data_mat <- generate_data(mechanism = mechanism, n = L, H = H, seed = seed + d)
    
    # Loop over the different training set dimensions 
    
    for (j in 1:length(Ns)){
      
      n_train <- Ns[j]
      M <- L - n_train
      
      train_mat <- data_mat[1:n_train,]
      train_list <- convert_features_list(train_mat)
      
      test_mat <- data_mat[(n_train+1):L, ]
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
      
      # Fill the expected value richness structures
      df_ev_richness_PoissonBB[d,j] <- mean(total_richness(object = PoissonBB_fit))
      df_ev_richness_NegBinBB[d,j] <- mean(total_richness(object = NegBinBB_fit))
      
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
      
      # A) Chao's rarefaction and extrapolation
      
      # Determine the frequency vector of the training sets
      Q_vec <- colSums(train_mat)
      Q_vec <- Q_vec[Q_vec>0]
      
      # Compute the curves with confidence intervals
      fit_Chao <- iNEXT.Sam(Spec = Q_vec, T = n_train, endpoint = L)
      
      rare_Chao <- as_tibble(fit_Chao[["q=0"]]) %>%
        select(-Cov.hat) %>%
        rename(medians = D0.hat, lbs = Norm.CI.Low, ubs = Norm.CI.High)
      tmp <- as.data.frame(rare_Chao)
      
      df_extr_last_Chao[d ,j ] <- tmp[nrow(tmp), "medians"]
      
      # B) Smoothed Good-Toulmin extrapolation
      sfs <- tabulate(colSums(train_mat))
      cts <- sapply(2:n_train, function(n) ncol(train_mat[1:n,colSums(train_mat[1:n,]) > 0])   )
      cts <- c(0, sum(train_mat[1,]) , cts)
      
      df_extr_last_GT[d, j] <- predict_good_toulmin(n_train, M, sfs, cts, alternative = 0)$preds [n_train+M+1]
      
      
    }
    
  }
  
  
  # Save the entire workspace related to the mechanism just performed
  save(list = ls(all.names = TRUE), file =  paste0("R_script_paper/",mechanism,"_fit_estimate_repeateddataset.RData"))
}



# 1) Main script single-dataset -------


## EFPF approach -------

# Choose mechanism
mechanism = "beta_pis" # c("homogeneous", "random_uniform", "broken_stick", "log-normal")

# Fit and estimate richness, rarefaction and extrapolation for GibbsFA's (save workspace)
if (!file.exists(paste0("R_script_paper/eb_EFPF_",mechanism,"_fit_estimate_singledataset.RData"))) {
  
  vars_fct_NegBinBB <- c(10, 1000) 
  vars_GammaIBP <- c(1, 1000)
  
  eb_init_BB <- list(alpha = -1, s = 100, Nhat_prime = 100)
  eb_known_BB <- list()
  
  eb_init_IBP <- list(alpha = 0.5, s = 1, Gamma = 10)
  eb_known_IBP <- list()
  
  # Call the routine to perform simulations
  eb_EFPF_fit_estimate_ecological_scenario_singledataset(mechanism = mechanism, 
                                                         vars_fct_NegBinBB,
                                                         vars_GammaIBP,
                                                         eb_init_BB, eb_known_BB,
                                                         eb_init_IBP, eb_known_IBP, 
                                                         seed = 123)
  
}

# Load the Work space
load(paste0("R_script_paper/eb_EFPF_",mechanism,"_fit_estimate_singledataset.RData"))

Ns
Kn <- sapply(Ns, function(n) sum(colSums(data_mat[1:n,]) > 0)  )
sum(colSums(data_mat) > 0)
plot(rarefaction(data_mat[1:1000,]))
print(paste0("features totali: ", H))
### Comparison Poisson, NegBin and Gamma ------

# 0.B) Check on rarefaction
n_rare <- Ns[2]
lab_comb_bb <- paste0("n_train.",n_rare,":Nbar.emp")
lab_comb_ibp <- paste0("n_train.",n_rare)

eb_EFPF_fit_PoissonBB_rare <- list_eb_EFPF_fit_PoissonBB[[lab_comb_bb]]
eb_EFPF_fit_NegBinBB_rare <- list_eb_EFPF_fit_NegBinBB[[1]][[lab_comb_bb]]
eb_EFPF_fit_GammaIBP_rare <- list_eb_EFPF_fit_GammaIBP[[1]][[lab_comb_ibp]]

accum_df <- tibble( x = 0:n_rare,
                    n_feat = c(0,rarefaction(data_mat[1:n_rare,], n_reorderings = 20)))

rare_EFPF_PoissonBB <- tibble( lambda_post = unname(unlist(
  rarefaction(object = eb_EFPF_fit_PoissonBB_rare, seed = seed)$lambda_post ))) %>%
  rename(means = lambda_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "PoissonBB",
             x = c(1:n_rare,0))


rare_EFPF_NegBinBB <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_EFPF_fit_NegBinBB_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "NegBinBB",
             x = c(1:n_rare,0))



rare_EFPF_GammaIBP <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_EFPF_fit_GammaIBP_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "GammaIBP",
             x = c(1:n_rare,0))

# for plot
rare_EFPF_mixtureBB <- rare_EFPF_NegBinBB %>%
  mutate(Model = "PoissonBB/NegBinBB")

df_rare <- rbind(rare_EFPF_mixtureBB,
                 rare_EFPF_GammaIBP)

df_rare$Model <- factor(df_rare$Model,
                        levels = c( "PoissonBB/NegBinBB", "GammaIBP"))

ggplot(df_rare, aes(x = x, y = means, color = Model)) +
  geom_line(linetype = "dashed", linewidth = 1) + 
  #facet_wrap(.~ Model, scales = "free_x", nrow = 1) +
  #geom_ribbon(aes(ymin = lb_bands, ymax = ub_bands), color = "red" , linewidth = 0.8, alpha = 0.1) +
  geom_point( data = accum_df, aes(x = x, y = n_feat), color="black", shape = 21, size = 0.05) +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) + 
  scale_color_tableau()
ggsave(filename = paste0("R_script_paper/Paper_plots/rarefaction_",mechanism,"_eb_EFPF.pdf"), width = 4, height = 4, dpi = 300, units = "in", device='pdf')


# 0.c) Check on K_n_r
n_knr <- Ns[2]
lab_comb_bb <- paste0("n_train.",n_knr,":Nbar.emp")
lab_comb_ibp <- paste0("n_train.",n_knr)

eb_EFPF_fit_PoissonBB_knr <- list_eb_EFPF_fit_PoissonBB[[lab_comb_bb]]
eb_EFPF_fit_NegBinBB_knr <- list_eb_EFPF_fit_NegBinBB[[1]][[lab_comb_bb]]
eb_EFPF_fit_GammaIBP_knr <- list_eb_EFPF_fit_GammaIBP[[1]][[lab_comb_ibp]]


observed_K_n_r <- tibble( r = 1:n_knr,
                          k_n_r = K_n_r(data_mat[1:n_knr,], n_reorderings = 1)[[paste0('N = ', n_knr)]])


K_n_r_EFPF_PoissonBB <- tibble( lambda_est = unname(unlist(
  K_n_r(object = eb_EFPF_fit_PoissonBB_knr, n = n_knr)[[paste0('N = ', n_knr)]]$lambda_est ))) %>%
  rename(means = lambda_est) %>%
  add_column(Model = "PoissonBB",
             r = 1:n_knr)


K_n_r_EFPF_NegBinBB <- tibble( mu0_est = unname(unlist(
  K_n_r(object = eb_EFPF_fit_NegBinBB_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
  rename(means = mu0_est) %>%
  add_column(Model = "NegBinBB",
             r = 1:n_knr)


K_n_r_EFPF_GammaIBP <- tibble( mu0_est = unname(unlist(
  K_n_r(object = eb_EFPF_fit_GammaIBP_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
  rename(means = mu0_est) %>%
  add_column(Model = "GammaIBP",
             r = 1:n_knr)

# for plot
K_n_r_EFPF_mixtureBB <- K_n_r_EFPF_PoissonBB %>%
  mutate(Model = "PoissonBB/NegBinBB")

df_K_n_r <- rbind(K_n_r_EFPF_mixtureBB,
                  K_n_r_EFPF_GammaIBP)

df_K_n_r$Model <- factor(df_K_n_r$Model, 
                         levels = c("PoissonBB/NegBinBB", "GammaIBP"))

r_positive <- observed_K_n_r %>%
  filter(k_n_r > 0) %>%
  select(r) %>%
  filter(r < 15)

df_K_n_r_plot <- df_K_n_r %>%
  filter(r %in% c(r_positive$r))

observed_K_n_r_plot <- observed_K_n_r %>%
  filter(r %in% c(r_positive$r))

ggplot(df_K_n_r_plot, aes(x = r, y = means, color = Model)) +
  geom_line( linetype = "dashed") +
  geom_point( data = observed_K_n_r_plot, aes(x = r, y = k_n_r), color="black", shape = 21, size = 1) +
  scale_y_log10() +
  xlab("r") + ylab(expression(m[r])) + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) + 
  scale_color_tableau()
ggsave(filename = paste0("R_script_paper/Paper_plots/knr_", mechanism, "_eb_EFPF.pdf"), width = 4, height = 4, dpi = 300, units = "in", device='pdf')



### Prediction -----------

#### Richness -------

# 1) Plot Richness: Point plot expected value (number of features)
labels_comb_bb <- paste(rep(paste("n_train", Ns, sep = "."), each = length(Nbars)+1),
                        c("Nbar.emp" , paste("Nbar", Nbars, sep = ".")), sep=":")

# PoissonBB
richness_EFPF_PoissonBB_df <- tibble(estimate = unname(sapply(list_eb_EFPF_fit_PoissonBB, function(x)
  total_richness(x)$lambda_post + ncol(x$feature_matrix)) ) ) %>%
  add_column(Model = "PoissonBB", 
             Nbar = rep(c("EB", Nbars), length(Ns)),
             n_train = rep(Ns, each = length(Nbars)+ 1)) 

# NegBinBB
richness_EFPF_NegBinBB_df <- tibble(estimate = numeric(), Model = character(),
                                    Nbar = character(), n_train = integer())

for (var_fct_NegBinBB in vars_fct_NegBinBB){
  list_eb_EFPF_fit_NegBinBB_var <- list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]]
  
  richness_EFPF_NegBinBB_df_var <- tibble(estimate = unname(sapply(list_eb_EFPF_fit_NegBinBB_var, function(x)
    total_richness(x)$mu0_post + ncol(x$feature_matrix)) ) ) %>%
    add_column(Model = paste0("NegBinBB x",var_fct_NegBinBB) , 
               Nbar = rep(c("EB", Nbars), length(Ns)),
               n_train = rep(Ns, each = length(Nbars)+ 1) )
  
  richness_EFPF_NegBinBB_df <- bind_rows(richness_EFPF_NegBinBB_df, richness_EFPF_NegBinBB_df_var)
  
}


joint_richness_long <- bind_rows(richness_EFPF_PoissonBB_df,
                                 richness_EFPF_NegBinBB_df) %>%
  mutate(Model = fct_relevel(Model, c("PoissonBB", 
                             paste0("NegBinBB x", vars_fct_NegBinBB))) )

n_train.labs <- paste0("n = ", Ns,", Kn = ", Kn )
names(n_train.labs) <- Ns

# plots estimator
ggplot(joint_richness_long, aes( y=estimate, x=Model, shape = Nbar)) +
  geom_point(size = 2) +
  facet_wrap(.~ n_train, labeller = labeller(n_train = n_train.labs ), scales = "free_x", nrow = 1) +
  theme_light() +
  geom_hline(aes(yintercept = H), linetype = "dashed") +
  theme(legend.position = "top") +
  ylab("Posterior mean of N") +
  scale_y_continuous(breaks = pretty_breaks()) +
  rremove("xlab") +
  scale_shape_discrete(name = "Prior mean of N") +
  theme(aspect.ratio = 1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = paste0("R_script_paper/Paper_plots/richness_points_", mechanism,"_eb_EFPF.pdf"), width = 9, height = 5, dpi = 300, units = "in", device='pdf')


# 2) Plot Richness: whole distributions for Nbar = 600 
Nbar_plot <- 400

# compute parametes of richness

# PoissonBB
list_eb_EFPF_est_PoissonBB <- list_eb_EFPF_fit_PoissonBB[grepl("Nbar.emp", names(list_eb_EFPF_fit_PoissonBB) )]
params_richness_EFPF_est_PoissonBB <- tibble( lambda_prime = unname(sapply(list_eb_EFPF_est_PoissonBB, function(x)
  total_richness(x)$lambda_post))) %>%
  add_column(n_train = Ns, Model = "PoissonBB", PM = "est") %>%
  mutate(lb = qpois(0.025, lambda_prime, lower.tail = TRUE, log.p = FALSE),
         ub = qpois(0.975, lambda_prime, lower.tail = TRUE, log.p = FALSE) )

list_eb_EFPF_fixed_PoissonBB <- list_eb_EFPF_fit_PoissonBB[grepl(paste0("Nbar.",Nbar_plot), names(list_eb_EFPF_fit_PoissonBB) )]
params_richness_EFPF_fixed_PoissonBB <- tibble( lambda_prime = unname(sapply(list_eb_EFPF_fixed_PoissonBB, function(x)
  total_richness(x)$lambda_post))) %>%
  add_column(n_train = Ns, Model = "PoissonBB", PM = "fixed") %>%
  mutate(lb = qpois(0.025, lambda_prime, lower.tail = TRUE, log.p = FALSE),
         ub = qpois(0.975, lambda_prime, lower.tail = TRUE, log.p = FALSE) )

# NegBinBB
params_richness_EFPF_est_NegBinBB <- tibble( n0_prime = numeric(),
                                             mu0_prime = numeric(),
                                             n_train = integer(),
                                             Model = character(),
                                             PM = character(),
                                             p_prime = numeric(),
                                             lb = numeric(), ub = numeric())
params_richness_EFPF_fixed_NegBinBB <- tibble( n0_prime = numeric(),
                                             mu0_prime = numeric(),
                                             n_train = integer(),
                                             Model = character(),
                                             PM = character(),
                                             p_prime = numeric(),
                                             lb = numeric(), ub = numeric())

for (var_fct_NegBinBB in vars_fct_NegBinBB){
  
  # est
  list_eb_EFPF_est_NegBinBB_var <- 
    list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]][grepl("Nbar.emp", names(list_eb_EFPF_fit_NegBinBB[[1]]) )]
  
  params_richness_EFPF_est_NegBinBB_var <- tibble( n0_prime = unname(sapply(list_eb_EFPF_est_NegBinBB_var, function(x)
    total_richness(x)$n0_post)),
    mu0_prime = unname(sapply(list_eb_EFPF_est_NegBinBB_var, function(x)
      total_richness(x)$mu0_post))) %>%
    add_column(n_train = Ns, Model = paste0("NegBinBB x", var_fct_NegBinBB),
               PM = "est") %>%
    mutate(p_prime = 1/(mu0_prime/n0_prime + 1),
           lb = qnbinom(0.025, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE),
           ub = qnbinom(0.975, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE) )
  
  params_richness_EFPF_est_NegBinBB <- bind_rows(params_richness_EFPF_est_NegBinBB,
                                                 params_richness_EFPF_est_NegBinBB_var)

  
  #fixed
  list_eb_EFPF_fixed_NegBinBB_var <- 
    list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]][grepl(paste0("Nbar.",Nbar_plot), names(list_eb_EFPF_fit_NegBinBB[[1]]) )]
  
  params_richness_EFPF_fixed_NegBinBB_var <- tibble( n0_prime = unname(sapply(list_eb_EFPF_fixed_NegBinBB_var, function(x)
    total_richness(x)$n0_post)),
    mu0_prime = unname(sapply(list_eb_EFPF_fixed_NegBinBB_var, function(x)
      total_richness(x)$mu0_post))) %>%
    add_column(n_train = Ns, Model = paste0("NegBinBB x", var_fct_NegBinBB),
               PM = "fixed") %>%
    mutate(p_prime = 1/(mu0_prime/n0_prime + 1),
           lb = qnbinom(0.025, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE),
           ub = qnbinom(0.975, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE) )
  
  
  params_richness_EFPF_fixed_NegBinBB <- bind_rows(params_richness_EFPF_fixed_NegBinBB,
                                                   params_richness_EFPF_fixed_NegBinBB_var)

  
  
}



bounds_distr <- tibble( lb = min(params_richness_EFPF_est_PoissonBB$lb + Kn, 
                                 params_richness_EFPF_est_NegBinBB$lb + Kn,
                                 params_richness_EFPF_fixed_PoissonBB$lb + Kn, 
                                 params_richness_EFPF_fixed_NegBinBB$lb + Kn,
                                 H),
                        ub = max(params_richness_EFPF_est_PoissonBB$ub + Kn,
                                 params_richness_EFPF_est_NegBinBB$ub + Kn,
                                 params_richness_EFPF_fixed_PoissonBB$ub + Kn,
                                 params_richness_EFPF_fixed_NegBinBB$ub + Kn,
                                 H))

# density computation
# PoissonBB
dens_richness_EFPF_est_PoissonBB <- bind_rows(lapply(1:length(Ns), function(j) 
  tibble( x = bounds_distr$lb: bounds_distr$ub) %>%
    mutate( y = dpois(x - Kn[j] , lambda = params_richness_EFPF_est_PoissonBB$lambda_prime[j])) %>%
    add_column(Model = "PoissonBB", PM = "est",
               n_train = params_richness_EFPF_est_PoissonBB$n_train[j])
))

dens_richness_EFPF_fixed_PoissonBB <- bind_rows(lapply(1:length(Ns), function(j) 
  tibble( x = bounds_distr$lb: bounds_distr$ub) %>%
    mutate( y = dpois(x - Kn[j] , lambda = params_richness_EFPF_fixed_PoissonBB$lambda_prime[j])) %>%
    add_column(Model = "PoissonBB", PM = "fixed",
               n_train = params_richness_EFPF_fixed_PoissonBB$n_train[j])
))

# NegBinBB
dens_richness_EFPF_est_NegBinBB <- tibble(x = integer(), y = numeric(),
                                          Model = character(),
                                          PM = character(),
                                          n_train = integer())

dens_richness_EFPF_fixed_NegBinBB <- tibble(x = integer(), y = numeric(),
                                            Model = character(),
                                            PM = character(),
                                            n_train = integer())


for (var_fct_NegBinBB in vars_fct_NegBinBB){
  
  # est
  params_richness_EFPF_est_NegBinBB_var <- params_richness_EFPF_est_NegBinBB %>%
    filter(Model == paste0("NegBinBB x", var_fct_NegBinBB))
  
  dens_richness_EFPF_est_NegBinBB_var <- bind_rows(lapply(1:length(Ns), function(j) 
    tibble( x = bounds_distr$lb: bounds_distr$ub) %>%
      mutate( y = dnbinom(x - Kn[j], size = params_richness_EFPF_est_NegBinBB_var$n0_prime[j], 
                          prob = params_richness_EFPF_est_NegBinBB_var$p_prime[j])) %>%
      add_column(Model = paste0("NegBinBB x", var_fct_NegBinBB), 
                 PM = "est",
                 n_train = params_richness_EFPF_est_NegBinBB_var$n_train[j])
  ))
  
  dens_richness_EFPF_est_NegBinBB <- bind_rows(dens_richness_EFPF_est_NegBinBB,
                                               dens_richness_EFPF_est_NegBinBB_var)
  
  
  # fixed
  params_richness_EFPF_fixed_NegBinBB_var <- params_richness_EFPF_fixed_NegBinBB %>%
    filter(Model == paste0("NegBinBB x", var_fct_NegBinBB))
  
  dens_richness_EFPF_fixed_NegBinBB_var <- bind_rows(lapply(1:length(Ns), function(j) 
    tibble( x = bounds_distr$lb: bounds_distr$ub) %>%
      mutate( y = dnbinom(x - Kn[j], size = params_richness_EFPF_fixed_NegBinBB_var$n0_prime[j], 
                          prob = params_richness_EFPF_fixed_NegBinBB_var$p_prime[j])) %>%
      add_column(Model = paste0("NegBinBB x", var_fct_NegBinBB),
                 PM = "fixed",
                 n_train = params_richness_EFPF_fixed_NegBinBB_var$n_train[j])
  ))
  
  dens_richness_EFPF_fixed_NegBinBB <- bind_rows(dens_richness_EFPF_fixed_NegBinBB,
                                                 dens_richness_EFPF_fixed_NegBinBB_var)
  
}


# Set Ns and Kn to plot
Ns_plot <- Ns
Kn_plot <- Kn

dens_richnesses <- rbind(dens_richness_EFPF_est_PoissonBB,
                         dens_richness_EFPF_est_NegBinBB,
                         dens_richness_EFPF_fixed_PoissonBB,
                         dens_richness_EFPF_fixed_NegBinBB ) %>%
  filter(n_train %in% Ns_plot)


dens_richnesses$Model <- factor(dens_richnesses$Model, 
                                levels = c("PoissonBB",
                                           paste0("NegBinBB x", vars_fct_NegBinBB)))

n_train.labs <- paste0("n = ", Ns_plot,", Kn = ", Kn_plot )
names(n_train.labs) <- Ns_plot
PM.labs <- c("est" = "EB", "fixed" = Nbar_plot )

ggplot(dens_richnesses, aes(x = x, y = y, color = Model)) +
  geom_line() +
  geom_vline(aes(xintercept = H), linetype="dashed") +
  facet_grid(PM ~ n_train,  
             labeller = labeller(n_train = n_train.labs, PM = PM.labs),
             scales = "free") +
  theme_light() +
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  xlab("# distinct features") + rremove("ylab") +
  scale_color_tableau() +
  theme(aspect.ratio = 1)
ggsave(filename = paste0("R_script_paper/Paper_plots/richness_distr_", mechanism, "_eb_EFPF.pdf"), width = 6, height = 6, dpi = 300, units = "in", device='pdf')



#### Extrapolation -----

# GT and Chao estimators

if (!file.exists(paste0("R_script_paper/eb_Freq_",mechanism,"_fit_estimate_singledataset.RData"))) {
  
  M <- 500
  
  # List to store Chao's estimates
  list_rare_extr_Chao <- vector(mode="list")
  # List to store smoothed GT's estimates
  list_extr_GT <- vector(mode="list")
  
  for (j in 1:length(Ns)){
    
    n_train <- Ns[j]
    
    train_mat <- data_mat[1:n_train,]
    train_mat <- train_mat[, colSums(train_mat) > 0]
    K <- ncol(train_mat)
    
    lab_comb_bb <- paste0("n_train.",n_train,":Nbar.emp")
    lab_comb_ibp <- paste("n_train", n_train, sep = ".")
    
    # A) Chao's rarefaction and extrapolation
    
    if (j < 3){
      # Determine the frequency vector of the training sets
      Q_vec <- colSums(train_mat)
      Q_vec <- Q_vec[Q_vec>0]
      
      # Compute the curves with confidence intervals
      fit_Chao <- iNEXT.Sam(Spec = Q_vec, T = n_train, endpoint = n_train + M)
      
      rare_Chao <- as_tibble(fit_Chao[["q=0"]]) %>%
        select(-Cov.hat) %>%
        rename(medians = D0.hat, lbs = Norm.CI.Low, ubs = Norm.CI.High)
      
      list_rare_extr_Chao[[lab_comb_ibp]] <- as.data.frame(rare_Chao)
      
    }
    
    
    # B) Smoothed Good-Toulmin extrapolation
    
    # Compute SFS vector and CTS vector
    sfs <- tabulate(colSums(train_mat))
    cts <- sapply(2:n_train, function(n) ncol(train_mat[1:n,colSums(train_mat[1:n,]) > 0])   )
    cts <- c(0, sum(train_mat[1,]) , cts)
    
    list_extr_GT[[lab_comb_ibp]] <- predict_good_toulmin(n_train, M, sfs, cts, alternative = 0)$preds
    
    
  }
  
  # Save the entire workspace related to the mechanism just performed
  save(list = ls(all.names = TRUE), file =  paste0("R_script_paper/eb_Freq_",mechanism,"_fit_estimate_singledataset.RData"))
  
}

# Load the Work space
load(paste0("R_script_paper/eb_Freq_",mechanism,"_fit_estimate_singledataset.RData"))



# Extract accumulation curve of the observed sample (or average accumulation)

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
list_eb_EFPF_PoissonBB <- list_eb_EFPF_fit_PoissonBB[grepl("Nbar.emp", names(list_eb_EFPF_fit_PoissonBB) )]
extr_EFPF_PoissonBB_df <- tibble(lambda = unname(unlist(lapply(list_eb_EFPF_PoissonBB, function(x)
  extrapolation(object = x, M = M, seed = seed)$lambda_post ))),
  n_train = rep(Ns, each = M),
  Kn = rep(Kn, each = M)) %>%
  mutate(lb = qpois(0.025, lambda, lower.tail = TRUE, log.p = FALSE),
         ub = qpois(0.975, lambda, lower.tail = TRUE, log.p = FALSE)) %>%
  rename(means = lambda) %>%
  add_row(means = 0, lb = 0, ub = 0, n_train = Ns, Kn = Kn) %>%
  mutate(means = means + Kn, lb = lb + Kn, ub = ub + Kn) %>%
  add_column(x = c(unlist(sapply(Ns, function(n) (n+1):(M+n))), Ns),
             Model = "PoissonBB")

extr_EFPF_PoissonBB_df$x <- as.integer(extr_EFPF_PoissonBB_df$x)
extr_EFPF_PoissonBB_df <- extr_EFPF_PoissonBB_df %>%
  select(means, lb, ub, n_train, x, Model)


# NegBin
extr_EFPF_NegBinBB_df <- tibble(means = numeric(), 
                                lb = numeric(), ub = numeric(),
                                n_train = integer(), 
                                x = integer(), Model = character())

for (var_fct_NegBinBB in vars_fct_NegBinBB){
  
  list_eb_EFPF_NegBinBB_var <- 
    list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]][grepl("Nbar.emp", names(list_eb_EFPF_fit_NegBinBB[[1]]) )]
  
  extr_EFPF_NegBinBB_df_var <- tibble(mu0 = unname(unlist(lapply(list_eb_EFPF_NegBinBB_var, function(x)
    extrapolation(object = x, M = M, seed = seed)$mu0_post ))),
    n0 = unname(unlist(lapply(list_eb_EFPF_NegBinBB_var, function(x)
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
               Model = paste0("NegBinBB x", var_fct_NegBinBB))
  
  extr_EFPF_NegBinBB_df_var$x <- as.integer(extr_EFPF_NegBinBB_df_var$x)
  extr_EFPF_NegBinBB_df_var <- extr_EFPF_NegBinBB_df_var %>%
    select(means, lb, ub, n_train, x, Model)
  
  extr_EFPF_NegBinBB_df <- bind_rows(extr_EFPF_NegBinBB_df, 
                                     extr_EFPF_NegBinBB_df_var)
  
}


# # GammaIBP
# extr_EFPF_GammaIBP_df <- tibble(means = numeric(), 
#                                 lb = numeric(), ub = numeric(),
#                                 n_train = integer(), 
#                                 x = integer(), Model = character())
# 
# for (var_GammaIBP in vars_GammaIBP){
#   
#   list_eb_EFPF_GammaIBP_var <- 
#     list_eb_EFPF_fit_GammaIBP[[paste0("var.", var_GammaIBP)]]
#   
#   extr_EFPF_GammaIBP_df_var <- tibble(mu0 = unname(unlist(lapply(list_eb_EFPF_GammaIBP_var, function(x)
#     extrapolation(object = x, M = M, seed = seed)$mu0_post ))),
#     n0 = unname(unlist(lapply(list_eb_EFPF_GammaIBP_var, function(x)
#       extrapolation(object = x, M = M, seed = seed)$n0_post ))),
#     n_train = rep(Ns, each = M),
#     Kn = rep(Kn, each = M)) %>%
#     mutate(p = 1/(mu0/n0 + 1),
#            lb = qnbinom(0.025, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE),
#            ub = qnbinom(0.975, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE)) %>%
#     rename(means = mu0) %>%
#     add_row(means = 0, lb = 0, ub = 0, n_train = Ns, Kn = Kn) %>%
#     mutate(means = means + Kn, lb = lb + Kn, ub = ub + Kn) %>%
#     add_column(x = c(unlist(sapply(Ns, function(n) (n+1):(M+n))), Ns),
#                Model = paste0("GammaIBP, var = ", var_GammaIBP))
#   
#   extr_EFPF_GammaIBP_df_var$x <- as.integer(extr_EFPF_GammaIBP_df_var$x)
#   extr_EFPF_GammaIBP_df_var <- extr_EFPF_GammaIBP_df_var %>%
#     select(means, lb, ub, n_train, x, Model)
#   
#   
#   extr_EFPF_GammaIBP_df <- bind_rows(extr_EFPF_GammaIBP_df, 
#                                      extr_EFPF_GammaIBP_df_var)
# 
# }

# Set Ns and Kn to plot
Ns_plot <- Ns[1:2]
Kn_plot <- Kn[1:2]

df_extr_GT_long <- list_extr_competitor_to_long(list_extr_GT, model = "GT") %>%
  rename(Model = model) %>%
  filter(t < n_train + M + 1, n_train %in% Ns_plot)
df_extr_Chao_long <- list_extr_competitor_to_long(list_rare_extr_Chao, model = "Chao") %>%
  rename(Model = model) %>%
  filter(t < n_train + M + 1, n_train %in% Ns_plot)

# Plot
temp <- tibble(n_train = Ns_plot, xvalues = Ns_plot)
extr_EFPF_PoissonBB_df <- extr_EFPF_PoissonBB_df %>%
  add_column(Model_gen = "PoissonBB/NegBinBB")
extr_EFPF_NegBinBB_df <- extr_EFPF_NegBinBB_df %>%
  add_column(Model_gen = "PoissonBB/NegBinBB")
# extr_EFPF_GammaIBP_df <- extr_EFPF_GammaIBP_df %>%
#   add_column(Model_gen = "GammaIBP")

extr_all_df <- rbind(extr_EFPF_PoissonBB_df, 
                     extr_EFPF_NegBinBB_df) %>%
  filter(n_train %in% Ns_plot)
                     

extr_all_df$Model <- factor(extr_all_df$Model, 
                            levels = c("PoissonBB", 
                                       paste0("NegBinBB x", vars_fct_NegBinBB),
                                       paste0("GammaIBP, var = ", vars_GammaIBP)))

accum_df_train <- accum_df_train %>%
  filter(n_train %in% Ns_plot)
accum_df_test <- accum_df_test %>%
  filter(n_train %in% Ns_plot)

n_train.labs <- paste0("n = ", Ns_plot,", Kn = ", Kn_plot )
names(n_train.labs) <- Ns_plot

ggplot(extr_all_df, aes(x, means, color = Model)) +
  geom_line(linetype = "dashed", linewidth = 1) +
  geom_ribbon(aes(ymin = lb, ymax = ub, color = Model), linewidth = 0.8, alpha = 0) +
  geom_point( data = accum_df_train, aes(x = x, y = n_feat),
              color="black", shape = 21, size = 1) +
  geom_point( data = accum_df_test, aes(x = x, y = n_feat),
              color="black", shape = 21, size = 0.5) +
  geom_line(data = df_extr_GT_long, aes(t, value), linetype = "dashed", linewidth = 0.9) +
  geom_line(data = df_extr_Chao_long, aes(t, medians), linetype = "dashed", linewidth = 0.9) +
  facet_grid(.~ n_train,
             labeller = labeller(n_train = n_train.labs),
             scales = "free_x")  +
  geom_vline(data = temp, mapping =  aes(xintercept = xvalues) , linetype = "dashed", color = "grey") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) +
  scale_color_tableau()
ggsave(filename = paste0("R_script_paper/Paper_plots/extr_", mechanism, "_eb_EFPF.pdf"), width = 9, height = 5, dpi = 300, units = "in", device='pdf')



# 2) Main script repeated-dataset -------


mechanism = "broken_stick" # c("homogeneous", "random_uniform", "broken_stick", "log-normal")
n_dataset = 2

# Fit and estimate richness, rarefaction and extrapolation for GibbsFA's (save workspace)
if (!file.exists(paste0("R_script_paper/",mechanism,"_fit_estimate_repeateddataset.RData"))) {
  
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
  
  fit_estimate_ecological_scenario_repeateddataset(mechanism = mechanism, n_dataset = n_dataset,
                                                   init_obj_PoissonBB, mcmcparams_obj_PoissonBB, prior_obj_PoissonBB,
                                                   init_obj_NegBinBB, mcmcparams_obj_NegBinBB, prior_obj_NegBinBB, c_fr,
                                                   init_obj_GammaIBP, mcmcparams_obj_GammaIBP, prior_obj_GammaIBP)
  
}

# Load the Work space
load(paste0("R_script_paper/",mechanism,"_fit_estimate_repeateddataset.RData"))


# 1) Plot richness boxplot 

df_ev_richness_PoissonBB_long <- df_ev_richness_PoissonBB %>% 
  pivot_longer(everything(), names_to = "training", values_to = "estimate") %>%
  add_column(model = "Poisson", n_train = rep(Ns[1:length(Ns)], each = n_dataset))

df_ev_richness_NegBinBB_long <- df_ev_richness_NegBinBB %>%
  pivot_longer(everything(), names_to = "training", values_to = "estimate") %>%
  add_column(model = "NegBin", n_train = rep(Ns[1:length(Ns)], each = n_dataset))

joint_richness_long <- bind_rows(df_ev_richness_PoissonBB_long, df_ev_richness_NegBinBB_long)

# plots
ggplot(joint_richness_long, aes(x=model, y=estimate)) +
  geom_boxplot() +
  facet_wrap(.~n_train, scales = "free_x", nrow = 1) +
  geom_hline(aes(yintercept = 500), linetype = "dashed") +
  theme_light() + 
  rremove("xlab") +
  ylab("# distinct features") +
  scale_y_continuous(breaks = pretty_breaks())+
  theme(aspect.ratio = 1)


# 2) Plot boxplots on accuracy 

df_new_PoissonBB <- df_extr_last_PoissonBB - df_obs_train
df_new_NegBinBB <- df_extr_last_NegBinBB - df_obs_train
df_new_GammaIBP <- df_extr_last_GammaIBP - df_obs_train
df_new_GT <- df_extr_last_GT - df_obs_train
df_new_Chao <- df_extr_last_Chao - df_obs_train

# compute accuracy
acc_alt_PoissonBB <- compute_accuracy(df_obs_new, df_new_PoissonBB, df_obs_train)
acc_alt_NegBinBB <- compute_accuracy(df_obs_new, df_new_NegBinBB, df_obs_train)
acc_alt_GammaIBP <- compute_accuracy(df_obs_new, df_new_GammaIBP, df_obs_train)
acc_alt_GT <- compute_accuracy(df_obs_new, df_new_GT, df_obs_train)
acc_alt_Chao <- compute_accuracy(df_obs_new, df_new_Chao, df_obs_train)

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

acc_alt_Chao_long <- acc_alt_Chao %>%
  pivot_longer( everything(), names_to = "training", values_to = "Accuracy") %>%
  add_column(model = "Chao", n_train = rep(Ns, each = n_dataset))


joint_alt_long <- bind_rows(acc_alt_PoissonBB_long, acc_alt_NegBinBB_long,
                            acc_alt_GammaIBP_long, acc_alt_GT_long, acc_alt_Chao_long)

# plots
ggplot(joint_alt_long, aes(x = model, y=Accuracy)) +
  geom_boxplot() +
  facet_wrap(.~n_train, scales = "free_x", nrow = 1) +
  theme_light() +
  rremove("xlab") +
  ylab("Error index") +
  scale_y_continuous(breaks = pretty_breaks())+
  theme(aspect.ratio = 1)


