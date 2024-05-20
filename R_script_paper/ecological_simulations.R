
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
    
    data_mat <- cbind(
      matrix(rbinom(n * 20, size = 1, prob = 0.9), n, 20),
      matrix(rbinom(n * 10, size = 1, prob = 0.8), n, 10),
      matrix(rbinom(n * 50, size = 1, prob = 0.7), n, 50),
      matrix(rbinom(n * 30, size = 1, prob = 0.6), n, 30),
      matrix(rbinom(n * 20, size = 1, prob = 0.5), n, 20),
      matrix(rbinom(n * 10, size = 1, prob = 0.4), n, 10),
      matrix(rbinom(n * 20, size = 1, prob = 0.3), n, 20),
      matrix(rbinom(n * 20, size = 1, prob = 0.2), n, 20),
      matrix(rbinom(n * 30, size = 1, prob = 0.1), n, 30),
      matrix(rbinom(n * 40, size = 1, prob = 0.05), n, 40),
      matrix(rbinom(n * 5, size = 1, prob = 0.01), n, 5))
    
    hist(colMeans(data_mat))
    
  } else if (mechanism == "beta_pis"){
    
    pis <- rbeta(H, 1, 20) # 1,20 is good
    
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


# EB Single dataset: function to produce fit and estimate on specific ecological scenario ------

eb_EFPF_fit_estimate_ecological_scenario_singledataset <- function(mechanism, 
                                                              eb_init_PoissonBB, eb_known_PoissonBB,
                                                              eb_init_NegBinBB, eb_known_NegBinBB, 
                                                              eb_init_GammaIBP, eb_known_GammaIBP,
                                                              seed = 1234){
  
  if (!mechanism %in% c("custom", "beta_pis", "homogeneous", "random_uniform", "broken_stick", "log-normal")){
    stop("Invalid generating mechanism.")
  }
  
  # Set training dimensions and dimension of the whole sample 
  Ns <- c(10000) #c(20, 40, 80) 
  N_max <- Ns[length(Ns)]
  
  if (mechanism == "log-normal"){
    L <- 300
  } else {
    L <- 10000
  }
  # Set maximum number of features
  H <- 500
  
  # Generate data 
  data_mat <- generate_data(mechanism = mechanism, n = L, H = H, seed = seed)
  
  # Set grid of desired value of E[N] = Nbar, with the empirical case to be added
  #Nbars <- c(300, 400, 600)
  
  # Structures to save fit and estimates 
  
  # List to store the fitted models for PoissonBB, NegBinBB and GammaIBP, under different settings 
  list_eb_fit_PoissonBB <- vector(mode = "list")
  list_eb_fit_NegBinBB <-  vector(mode = "list")
  list_eb_fit_GammaIBP <-  vector(mode = "list")
  
  # eb_NegBinBB_list_variances <- vector(mode="list", length = length(cfrs))
  # names(eb_NegBinBB_list_variances) <- paste0("c_fr.", cfrs)
   
  # List of Nbar_emp for the different training sets
  list_Nbar_emp <- vector(mode = "list")
  
  # List to store Chao's estimates
  list_rare_extr_Chao <- vector(mode="list")
  # List to store smoothed GT's estimates
  list_extr_GT <- vector(mode="list")
  
  # for (cfr in cfrs){
  #   eb_known_NegBinBB_vars <- list(var_fct = cfr)
  #   eb_params_obj_NegBinBB <- eb_params(model = "NegBinBB",
  #                                       init = eb_init_NegBinBB,
  #                                       known = eb_known_NegBinBB_vars)
  #   
  #   eb_NegBinBB_list[[paste0("c_fr.", cfr)]] <- GibbsFA_eb(feature_matrix = data_mat[1:N_max,],
  #                                                          model = "NegBinBB",
  #                                                          eb_params =  eb_params_obj_NegBinBB)
  # }
  
  # Loop over the different training set dimensions 
  
  for (j in 1:length(Ns)){
    
    n_train <- Ns[j]
    M <- L - n_train
    
    train_mat <- data_mat[1:n_train,]
    
    lab_comb_bb <- paste0("n_train.",n_train,":Nbar.emp")
    lab_comb_ibp <- paste("n_train", n_train, sep = ".")
    
    
    # Fit the models
    # PoissonBB
    eb_params_obj_PoissonBB <- eb_params(model = "PoissonBB", 
                                         init = eb_init_PoissonBB, 
                                         known = eb_known_PoissonBB)
    
    eb_PoissonBB_fit <- GibbsFA_eb(feature_matrix = train_mat, 
                                   model = "PoissonBB", type = "EFPF",
                                   eb_params =  eb_params_obj_PoissonBB)
    
    # NegBinBB
    eb_params_obj_NegBinBB <- eb_params(model = "NegBinBB",
                                        init = eb_init_NegBinBB,
                                        known = eb_known_NegBinBB)

    eb_NegBinBB_fit <- GibbsFA_eb(feature_matrix = train_mat,
                                  model = "NegBinBB", type = "EFPF",
                                  eb_params =  eb_params_obj_NegBinBB)

    # # GammaIBP
    # eb_params_obj_GammaIBP <- eb_params(model = "GammaIBP",
    #                                     init = eb_init_GammaIBP,
    #                                     known = eb_known_GammaIBP)
    # 
    # eb_GammaIBP_fit <- GibbsFA_eb(feature_matrix = train_mat,
    #                               model = "GammaIBP",
    #                               eb_params =  eb_params_obj_GammaIBP)

    # Fill the structures of MCMC chains of the parameters
    list_eb_fit_PoissonBB[[lab_comb_bb]] <- eb_PoissonBB_fit
    list_eb_fit_NegBinBB[[lab_comb_bb]] <- eb_NegBinBB_fit
    #list_eb_fit_GammaIBP[[lab_comb_bb]] <- eb_GammaIBP_fit


    # Loop over the different Nbar hyperparameters (only for BB mixtures) 
    
    # for (v in 1:length(Nbars)){
    # 
    #   Nbar <- Nbars[v]
    # 
    #   # Fit the models
    #   # PoissonBB
    #   eb_known_PoissonBB_Nbar <- list("lambda" = Nbar)
    #   eb_init_PoissonBB_Nbar <- eb_init_PoissonBB[! names(eb_init_PoissonBB) %in% c("lambda")]
    #   eb_params_obj_PoissonBB <- eb_params(model = "PoissonBB",
    #                                        init = eb_init_PoissonBB_Nbar,
    #                                        known = eb_known_PoissonBB_Nbar)
    # 
    #   eb_PoissonBB_fit <- GibbsFA_eb(feature_matrix = train_mat,
    #                                  model = "PoissonBB",
    #                                  eb_params =  eb_params_obj_PoissonBB)
    # 
    #   # NegBinBB
    #   eb_known_NegBinBB_Nbar <- list("var_fct" = c_fr,
    #                                  "mu0" = Nbar)
    #   eb_init_NegBinBB_Nbar <- eb_init_NegBinBB[! names(eb_init_NegBinBB) %in% c("mu0", "var_fct")]
    #   eb_params_obj_NegBinBB <- eb_params(model = "NegBinBB",
    #                                       init = eb_init_NegBinBB_Nbar,
    #                                       known = eb_known_NegBinBB_Nbar)
    # 
    #   eb_NegBinBB_fit <- GibbsFA_eb(feature_matrix = train_mat,
    #                                 model = "NegBinBB",
    #                                 eb_params =  eb_params_obj_NegBinBB)
    # 
    #   # Fill the structures of MCMC chains of the parameters
    #   lab_comb_bb <- paste0("n_train.",n_train,":Nbar.", Nbar)
    # 
    #   list_eb_fit_PoissonBB[[lab_comb_bb]] <- eb_PoissonBB_fit
    #   list_eb_fit_NegBinBB[[lab_comb_bb]] <- eb_NegBinBB_fit
    # 
    # 
    # }
    # 
    # # Competitors
    # 
    # # A) Chao's rarefaction and extrapolation
    # 
    # # Determine the frequency vector of the training sets
    # Q_vec <- colSums(train_mat)
    # Q_vec <- Q_vec[Q_vec>0]
    # 
    # # Compute the curves with confidence intervals
    # fit_Chao <- iNEXT.Sam(Spec = Q_vec, T = n_train, endpoint = L)
    # 
    # rare_Chao <- as_tibble(fit_Chao[["q=0"]]) %>%
    #   select(-Cov.hat) %>%
    #   rename(medians = D0.hat, lbs = Norm.CI.Low, ubs = Norm.CI.High)
    # 
    # list_rare_extr_Chao[[lab_comb_ibp]] <- as.data.frame(rare_Chao)
    # 
    # 
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
  save(list = ls(all.names = TRUE), file =  paste0("R_script_paper/eb_EFPF_",mechanism,"_fit_estimate_singledataset.RData"))
}



 eb_MM_fit_estimate_ecological_scenario_singledataset <- function(mechanism, 
                                                                 var_fct, #cfrs,
                                                                 eb_init_GammaIBP, eb_known_GammaIBP,
                                                                 seed = 1234){
  
  if (!mechanism %in% c("custom", "beta_pis", "homogeneous", "random_uniform", "broken_stick", "log-normal")){
    stop("Invalid generating mechanism.")
  }
  
  # Set training dimensions and dimension of the whole sample 
  Ns <- c(50, 500, 5000) 
  N_max <- Ns[length(Ns)]
  
  if (mechanism == "log-normal"){
    L <- 600
  } else {
    L <- N_max + 300
  }
  # Set maximum number of features
  H <- 500
  
  # Generate data 
  data_mat <- generate_data(mechanism = mechanism, n = L, H = H, seed = seed)
  
  # Set grid of desired value of E[N] = Nbar, with the empirical case to be added
  Nbars <- c(300, 400, 600)
  
  # Structures to save fit and estimates 
  
  # List to store the fitted models for PoissonBB, NegBinBB and GammaIBP, under different settings 
  list_eb_MM_biased_fit_PoissonBB <- vector(mode = "list")
  list_eb_MM_biased_fit_NegBinBB <-  vector(mode = "list")
  list_eb_MM_cens_fit_PoissonBB <- vector(mode = "list")
  list_eb_MM_cens_fit_NegBinBB <-  vector(mode = "list")
  list_eb_MM_fit_GammaIBP <-  vector(mode = "list")

  # Loop over the different training set dimensions 
  
  for (j in 1:length(Ns)){
    
    n_train <- Ns[j]
    M <- L - n_train
    
    train_mat <- data_mat[1:n_train,]
    
    lab_comb_bb <- paste0("n_train.",n_train,":Nbar.emp")
    lab_comb_ibp <- paste("n_train", n_train, sep = ".")
    
    # Fit the models
    # PoissonBB
    # eb_MM_biased_PoissonBB_fit <- GibbsFA_eb(feature_matrix = train_mat,
    #                                model = "PoissonBB", 
    #                                type = "MM_biased")
    # 
    eb_MM_cens_PoissonBB_fit <- GibbsFA_eb(feature_matrix = train_mat, 
                                           model = "PoissonBB", 
                                           type = "MM_censored")
    
    # NegBinBB
    # eb_MM_biased_NegBinBB_fit <- GibbsFA_eb(feature_matrix = train_mat,
    #                                          model = "NegBinBB", 
    #                                          type = "MM_biased", var_fct = var_fct)
    # 
    eb_MM_cens_NegBinBB_fit <- GibbsFA_eb(feature_matrix = train_mat, 
                                           model = "NegBinBB", 
                                           type = "MM_censored", var_fct = var_fct)
    
    # # GammaIBP
    # eb_params_obj_GammaIBP <- eb_params(model = "GammaIBP",
    #                                     init = eb_init_GammaIBP,
    #                                     known = eb_known_GammaIBP)
    # 
    # eb_MM_GammaIBP_fit <- GibbsFA_eb(feature_matrix = train_mat,
    #                                  model = "GammaIBP", type = "MM",
    #                                  eb_params =  eb_params_obj_GammaIBP)
    # 
    # Fill the structures of MCMC chains of the parameters
    #list_eb_MM_biased_fit_PoissonBB[[lab_comb_bb]] <- eb_MM_biased_PoissonBB_fit
    list_eb_MM_cens_fit_PoissonBB[[lab_comb_bb]] <- eb_MM_cens_PoissonBB_fit
    
    #list_eb_MM_biased_fit_NegBinBB[[lab_comb_bb]] <- eb_MM_biased_NegBinBB_fit
    list_eb_MM_cens_fit_NegBinBB[[lab_comb_bb]] <- eb_MM_cens_NegBinBB_fit
    
    #list_eb_MM_fit_GammaIBP[[lab_comb_bb]] <- eb_MM_GammaIBP_fit
    
    
    # Loop over the different Nbar hyperparameters (only for BB mixtures) 
    
    for (v in 1:length(Nbars)){

      Nbar <- Nbars[v]

      # Fit the models
      # PoissonBB
      # eb_MM_biased_PoissonBB_fit <- GibbsFA_eb(feature_matrix = train_mat,
      #                                          model = "PoissonBB",
      #                                          type = "MM_biased", Nhat_MM = Nbar)

      eb_MM_cens_PoissonBB_fit <- GibbsFA_eb(feature_matrix = train_mat,
                                             model = "PoissonBB",
                                             type = "MM_censored", Nhat_MM = Nbar)

      # NegBinBB
      # eb_MM_biased_NegBinBB_fit <- GibbsFA_eb(feature_matrix = train_mat,
      #                                         model = "NegBinBB",
      #                                         type = "MM_biased", Nhat_MM = Nbar,
      #                                         var_fct = var_fct)

      eb_MM_cens_NegBinBB_fit <- GibbsFA_eb(feature_matrix = train_mat,
                                            model = "NegBinBB",
                                            type = "MM_censored", Nhat_MM = Nbar,
                                            var_fct = var_fct)

      # Fill the structures of MCMC chains of the parameters
      lab_comb_bb <- paste0("n_train.",n_train,":Nbar.", Nbar)

      #list_eb_MM_biased_fit_PoissonBB[[lab_comb_bb]] <- eb_MM_biased_PoissonBB_fit
      list_eb_MM_cens_fit_PoissonBB[[lab_comb_bb]] <- eb_MM_cens_PoissonBB_fit

      #list_eb_MM_biased_fit_NegBinBB[[lab_comb_bb]] <- eb_MM_biased_NegBinBB_fit
      list_eb_MM_cens_fit_NegBinBB[[lab_comb_bb]] <- eb_MM_cens_NegBinBB_fit


    }

  }
  
  # Save the entire workspace related to the mechanism just performed
  save(list = ls(all.names = TRUE), file =  paste0("R_script_paper/eb_MM_",mechanism,"_fit_estimate_singledataset.RData"))
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



########### 1) Main script single-dataset -------

### 1.A) Prior approach ----

# Choose mechanism
mechanism = "broken_stick" # c("homogeneous", "random_uniform", "broken_stick", "log-normal")
  
# Fit and estimate richness, rarefaction and extrapolation for GibbsFA's (save workspace)
if (!file.exists(paste0("R_script_paper/",mechanism,"_fit_estimate_singledataset.RData"))) {
  
  # 1) PoissonBB 
  
  # Initialization and MCMC setting 
  init_PoissonBB <- list(alpha_0 = -1, s_0 = 1)
  init_obj_PoissonBB <- initialization(model = "PoissonBB", init = init_PoissonBB )
  mcmcparams_PoissonBB <- list(tau = 0.1, S = 3*10^4, n_burnin = 10^3, thin = 2)
  mcmcparams_obj_PoissonBB <- mcmcparameters(model = "PoissonBB", mcmcparams = mcmcparams_PoissonBB)
  
  # Hyperparameters elicitation 
  hyper_PoissonBB <- list(a_alpha = 0.001, b_alpha = 0.0001,
                          a_s = 0.001, b_s = 0.0001,
                          lambda = 1)
  prior_obj_PoissonBB <- prior(model = "PoissonBB", hyper = hyper_PoissonBB) 
  summary(prior_obj_PoissonBB)
  
  # 2) NegBinBB
  
  # Initialization and MCMC setting
  init_NegBinBB <- list(alpha_0 = -1, s_0 = 1)
  init_obj_NegBinBB <- initialization(model = "NegBinBB", init = init_NegBinBB )
  mcmcparams_NegBinBB <- list(tau = 0.1, S = 3*10^4, n_burnin = 10^3, thin = 2)
  mcmcparams_obj_NegBinBB <- mcmcparameters(model = "NegBinBB", mcmcparams = mcmcparams_NegBinBB)
  
  # Hyperparameters elicitation 
  c_fr <- 10
  
  hyper_NegBinBB <- list(a_alpha = 0.001, b_alpha = 0.0001,
                         a_s = 0.001, b_s = 0.0001,
                         n0 = 1, # n0, mu0 are set s.t. E(N) = Nbar, Var(N) = c_fr*E(N)
                         mu0 = 0.5)
  prior_obj_NegBinBB <- prior(model = "NegBinBB", hyper = hyper_NegBinBB) 
  summary(prior_obj_NegBinBB)
  
  # 3) GammaIBP
  
  # Initialization and MCMC setting 
  init_GammaIBP <- list(alpha_0 = 0.5, s_0 = 15, a_0 = 5, b_0 = 1)
  init_obj_GammaIBP <- initialization(model = "GammaIBP", init = init_GammaIBP )
  mcmcparams_GammaIBP <- list(sigq_alpha = 0.1, sigq_s = 0.1, 
                              S = 3*10^4, n_burnin = 10^3, thin = 2)
  mcmcparams_obj_GammaIBP <- mcmcparameters(model = "GammaIBP", mcmcparams = mcmcparams_GammaIBP)
  
  # Hyperparameters elicitation 
  hyper_GammaIBP <- list(a_alpha = 0.0001, b_alpha = 0.0001,
                         a_s = 0.001, b_s = 0.0001,
                         q = 0.05, r = 0.1, t = 0.01)
  prior_obj_GammaIBP <- prior(model = "GammaIBP", hyper = hyper_GammaIBP) 
  summary(prior_obj_GammaIBP)
  
  # 4) Call the routine to perform simulations
  fit_estimate_ecological_scenario_singledataset(mechanism = mechanism, 
                                                 init_obj_PoissonBB, mcmcparams_obj_PoissonBB, prior_obj_PoissonBB,
                                                 init_obj_NegBinBB, mcmcparams_obj_NegBinBB, prior_obj_NegBinBB, c_fr,
                                                 init_obj_GammaIBP, mcmcparams_obj_GammaIBP, prior_obj_GammaIBP )
  
}

# Load the Work space
load(paste0("R_script_paper/",mechanism,"_fit_estimate_singledataset.RData"))

# 1) Plot Richness: Point plot expected value (number of features)
labels_comb_bb <- paste(rep(paste("n_train", Ns, sep = "."), each = length(Nbars)+1),
                        c("Nbar.emp" , paste("Nbar", Nbars, sep = ".")), sep=":")

richness_PoissonBB_long <- gather(as_tibble(list_richness_PoissonBB), training, estimate,
                                  labels_comb_bb,
                                  factor_key=TRUE) %>%
  add_column(Model= "PoissonBB",
             n_train = rep(Ns, each = length(list_richness_PoissonBB[[1]])*(length(Nbars)+1)),
             Nbar = rep(rep(c("emp", Nbars), each = length(list_richness_PoissonBB[[1]])), length(Ns) ) )

richness_NegBinBB_long <- gather(as_tibble(list_richness_NegBinBB), training, estimate,
                                labels_comb_bb,
                                factor_key=TRUE) %>%
  add_column(Model= "NegBinBB",
             n_train = rep(Ns, each = length(list_richness_NegBinBB[[1]])*(length(Nbars)+1)),
             Nbar = rep(rep(c("emp", Nbars), each = length(list_richness_NegBinBB[[1]])), length(Ns) ))

joint_richness_long <- bind_rows(richness_PoissonBB_long, richness_NegBinBB_long) %>%
  mutate(Model = fct_relevel(Model, c( "NegBinBB", "PoissonBB")))

table_richness <- joint_richness_long %>% group_by(n_train, Nbar, Model) %>%
  summarise(estimator = mean(estimate),
            lb = quantile(estimate, probs = 0.025),
            ub = quantile(estimate, probs = 0.975), .groups = 'drop') %>%
  mutate(Model_spec = paste0(Model,":Nbar.",Nbar) )

# plots estimator
dev.new()
ggplot(table_richness, aes( y=estimator, x=Model, shape = Nbar)) +
  geom_point() +
  facet_wrap(.~n_train, scales = "free_x", nrow = 1) +
  theme_light() +
  geom_hline(aes(yintercept = 500), linetype = "dashed") +
  theme(legend.position = "top") +
  ylab("Posterior mean of N") +
  scale_y_continuous(breaks = pretty_breaks()) +
  rremove("xlab") +
  scale_shape_discrete(name = "Prior mean of N") +
  theme(aspect.ratio = 1)


# 2) Plot Richness: whole distributions for Nbar = 400
joint_400_long <- joint_richness_long %>%
  filter(Nbar == "400") %>%
  select(-training)

dev.new()
ggplot(joint_400_long, aes(x = estimate, color = Model)) +
  stat_density(aes(x=estimate, colour=Model),
               geom="line",position="identity", adjust = 2.5) +
  geom_vline(aes(xintercept = 500), linetype="dashed") +
  facet_wrap(.~"n = "n_train, scales = "free_x") +
  theme_light() +
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  xlab("# distinct features") + rremove("ylab") +
  scale_color_tableau() +
  theme(aspect.ratio = 1)


# 3) Plot Extrapolation

df_extr_PoissonBB_long <- list_extr_GibbsFA_to_long(list_extr_PoissonBB, model = "Poisson" )
df_extr_NegBinBB_long <- list_extr_GibbsFA_to_long(list_extr_NegBinBB, model = "NegBin")
df_extr_GammaIBP_long <- list_extr_GibbsFA_to_long(list_extr_GammaIBP, model = "Gamma")

df_extr_GT_long <- list_extr_competitor_to_long(list_extr_GT, model = "GT")  
df_extr_Chao_long <- list_extr_competitor_to_long(list_rare_extr_Chao, model = "Chao")  

accum <- rarefaction(data_mat)
accum_df <- data.frame("accum" = c(0, accum),
                       "t" = 0:length(accum))

df_GibbsFA <- rbind(df_extr_PoissonBB_long, df_extr_NegBinBB_long, df_extr_GammaIBP_long) %>%
  filter(Nbar %in% c("Not applicable", "emp"))

temp <- tibble(n_train = Ns, xvalues = Ns)

dev.new()
ggplot(df_GibbsFA, aes(x = t, y = means, color = model)) +
  geom_line() +
  geom_line(data = accum_df, aes(t, accum), color="black", linetype="solid") +
  geom_line(data = df_extr_GT_long, aes(t, value)) +
  geom_line(data = df_extr_Chao_long, aes(t, medians)) +
  facet_wrap(.~ n_train, scales = "free_x") +
  geom_vline(data = temp, mapping =  aes(xintercept = xvalues) , linetype = "dashed", color = "grey") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_tableau() +
  theme(aspect.ratio = 1)




### 1.B) EFPF approach ----

# Choose mechanism
mechanism = "beta_pis" # c("homogeneous", "random_uniform", "broken_stick", "log-normal")

# Fit and estimate richness, rarefaction and extrapolation for GibbsFA's (save workspace)
if (!file.exists(paste0("R_script_paper/eb_EFPF_",mechanism,"_fit_estimate_singledataset.RData"))) {
  
  # 1) PoissonBB 
  
  # Initialization and known parameters
  eb_init_PoissonBB <- list(alpha = -10, s = 1, lambda = 500)
  eb_known_PoissonBB <- list()
  
  # 2) NegBinBB
  
  # Initialization and known parameters
  c_fr <- 5
  eb_init_NegBinBB <- list(alpha = - 100, s = 100, mu0 = 500)
  eb_known_NegBinBB <- list(var_fct = c_fr)
  
  # 3) GammaIBP
  
  # Initialization and known parameters
  eb_init_GammaIBP <- list(alpha = 0.5, s = 1, a = 1, b = 1)
  eb_known_GammaIBP <- list()
  
  # NegBin for different variances
  #cfrs <- c(2, 10, 50)
  
  # 4) Call the routine to perform simulations
  eb_EFPF_fit_estimate_ecological_scenario_singledataset(mechanism = mechanism, 
                                                    eb_init_PoissonBB, eb_known_PoissonBB,
                                                    eb_init_NegBinBB, eb_known_NegBinBB,
                                                    eb_init_GammaIBP, eb_known_GammaIBP, seed = 1234)
  
}

# Load the Work space
load(paste0("R_script_paper/eb_EFPF_",mechanism,"_fit_estimate_singledataset.RData"))

Kn <- sapply(Ns, function(n) sum(colSums(data_mat[1:n,]) > 0)  )
sum(colSums(data_mat) > 0)

###### Comparison Poisson, NegBin and Gamma ------

# 0) Model checking: verify model is compatible with data 
plot( x = 0:100, y = c(0, rarefaction(data_mat)[1:100]) )

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



# 0.B) Check on rarefaction
n_rare <- Ns[3]
eb_fit_PoissonBB_rare <- list_eb_fit_PoissonBB[[2]]
eb_fit_NegBinBB_rare <- list_eb_fit_NegBinBB[[3]]
eb_fit_GammaIBP_rare <- list_eb_fit_GammaIBP[[3]]

accum_df <- tibble( x = 0:n_rare,
                    n_feat = c(0,rarefaction(data_mat[1:n_rare,], n_reorderings = 10)))

rare_PoissonBB <- tibble( lambda_post = unname(unlist(
  rarefaction(object = eb_fit_PoissonBB_rare, seed = seed)$lambda_post ))) %>%
  rename(means = lambda_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "PoissonBB",
             x = c(1:n_rare,0))

rare_NegBinBB <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_fit_NegBinBB_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "NegBinBB",
             x = c(1:n_rare,0))

rare_GammaIBP <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_fit_GammaIBP_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "GammaIBP",
             x = c(1:n_rare,0)) 

df_rare <- rbind(rare_PoissonBB, rare_NegBinBB, rare_GammaIBP)
df_rare$Model <- factor(df_rare$Model, levels = c("PoissonBB", "NegBinBB", "GammaIBP"))

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
ggsave(filename = "R_script_paper/Paper_plots/rarefaction_beta_pis_eb.pdf", width = 8, height = 3.5, dpi = 300, units = "in", device='pdf')


# 0.c) Check on K_n_r
n_knr <- Ns[3]
eb_fit_PoissonBB_knr <- list_eb_fit_PoissonBB[[9]]
eb_fit_NegBinBB_knr <- list_eb_fit_NegBinBB[[9]]
eb_fit_GammaIBP_knr <- list_eb_fit_GammaIBP[[3]]


observed_K_n_r <- tibble( r = 1:n_knr,
                          k_n_r = K_n_r(data_mat[1:n_knr,], n_reorderings = 1)[[paste0('N = ', n_knr)]])

K_n_r_PoissonBB <- tibble( lambda_est = unname(unlist(
  K_n_r(object = eb_fit_PoissonBB_knr, n = n_knr)[[paste0('N = ', n_knr)]]$lambda_est ))) %>%
  rename(means = lambda_est) %>%
  add_column(Model = "PoissonBB",
             r = 1:n_knr)

K_n_r_NegBinBB <- tibble( mu0_est = unname(unlist(
  K_n_r(object = eb_fit_NegBinBB_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
  rename(means = mu0_est) %>%
  add_column(Model = "NegBinBB",
             r = 1:n_knr)

K_n_r_GammaIBP <- tibble( mu0_est = unname(unlist(
  K_n_r(object = eb_fit_GammaIBP_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
  rename(means = mu0_est) %>%
  add_column(Model = "GammaIBP",
             r = 1:n_knr)

df_K_n_r <- rbind(K_n_r_PoissonBB, K_n_r_NegBinBB, K_n_r_GammaIBP)
df_K_n_r$Model <- factor(df_K_n_r$Model, levels = c("PoissonBB", "NegBinBB", "GammaIBP"))

r_positive <- observed_K_n_r %>%
  filter(k_n_r > 0) %>%
  select(r)

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
ggsave(filename = "R_script_paper/Paper_plots/knr_beta_pis_eb.pdf", width = 5, height = 4.5, dpi = 300, units = "in", device='pdf')


##### Comparison among NegBin with different variances -----------
# 0.B) Check on rarefaction
n_rare <- N_max
eb_fit_NegBinBB_1_rare <- eb_NegBinBB_list[[1]]
eb_fit_NegBinBB_2_rare <- eb_NegBinBB_list[[2]]
eb_fit_NegBinBB_3_rare <- eb_NegBinBB_list[[3]]

accum_df <- tibble( x = 0:n_rare,
                    n_feat = c(0,rarefaction(data_mat[1:n_rare,], n_reorderings = 10)))



rare_NegBinBB_1 <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_fit_NegBinBB_1_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Variance = paste0(cfrs[1]," x"),
             x = c(1:n_rare,0))

rare_NegBinBB_2 <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_fit_NegBinBB_2_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Variance = paste0(cfrs[2]," x"),
             x = c(1:n_rare,0))

rare_NegBinBB_3 <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_fit_NegBinBB_3_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Variance = paste0(cfrs[3]," x"),
             x = c(1:n_rare,0))

df_rare <- rbind(rare_NegBinBB_1, rare_NegBinBB_2, rare_NegBinBB_3)
df_rare$Variance <- factor(df_rare$Variance, 
                           levels = paste0(cfrs, " x"))

ggplot(df_rare, aes(x = x, y = means, color = Variance)) +
  geom_line(linetype = "solid", color = "red" , linewidth = 0.9) +
  facet_wrap(.~ Variance, scales = "free_x", nrow = 1) +
  #geom_ribbon(aes(ymin = lb_bands, ymax = ub_bands), color = "red" , linewidth = 0.8, alpha = 0.1) +
  geom_point( data = accum_df, aes(x = x, y = n_feat), color="black", shape = 1, size = 0.5) +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1)
ggsave(filename = "R_script_paper/Paper_plots/rarefaction__variances_beta_pis_eb.pdf", width = 8, height = 3.5, dpi = 300, units = "in", device='pdf')


# 0.c) Check on K_n_r
n_knr <- N_max
eb_fit_NegBinBB_1_knr <- eb_NegBinBB_list[[1]]
eb_fit_NegBinBB_2_knr <- eb_NegBinBB_list[[2]]
eb_fit_NegBinBB_3_knr <- eb_NegBinBB_list[[3]]


observed_K_n_r <- tibble( r = 1:n_knr,
                          k_n_r = K_n_r(data_mat[1:n_knr,], n_reorderings = 1)[[paste0('N = ', n_knr)]])

K_n_r_NegBinBB_1 <- tibble( mu0_est = unname(unlist(
  K_n_r(object = eb_fit_NegBinBB_1_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
  rename(means = mu0_est) %>%
  add_column(Variance = paste0(cfrs[1], " x"),
             r = 1:n_knr)

K_n_r_NegBinBB_2 <- tibble( mu0_est = unname(unlist(
  K_n_r(object = eb_fit_NegBinBB_2_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
  rename(means = mu0_est) %>%
  add_column(Variance = paste0(cfrs[2], " x"),
             r = 1:n_knr)

K_n_r_NegBinBB_3 <- tibble( mu0_est = unname(unlist(
  K_n_r(object = eb_fit_NegBinBB_3_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
  rename(means = mu0_est) %>%
  add_column(Variance = paste0(cfrs[3], " x"),
             r = 1:n_knr)

df_K_n_r <- rbind(K_n_r_NegBinBB_1, K_n_r_NegBinBB_2, K_n_r_NegBinBB_3)
df_K_n_r$Variance <- factor(df_K_n_r$Variance, 
                            levels = paste0(cfrs, " x"))

r_positive <- observed_K_n_r %>%
  filter(k_n_r > 0) %>%
  select(r)

df_K_n_r_plot <- df_K_n_r %>%
  filter(r %in% c(r_positive$r))

observed_K_n_r_plot <- observed_K_n_r %>%
  filter(r %in% c(r_positive$r))

ggplot(df_K_n_r_plot, aes(x = r, y = means, color = Variance)) +
  geom_line( linetype = "dashed") +
  geom_point( data = observed_K_n_r_plot, aes(x = r, y = k_n_r), color="black", size = 1) +
  scale_y_log10() +
  xlab("r") + ylab(expression(m[r])) + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1)
ggsave(filename = "R_script_paper/Paper_plots/knr_variances_beta_pis_eb.pdf", width = 5, height = 4.5, dpi = 300, units = "in", device='pdf')




# 1) Plot Richness: Point plot expected value (number of features)
labels_comb_bb <- paste(rep(paste("n_train", Ns, sep = "."), each = length(Nbars)+1),
                        c("Nbar.emp" , paste("Nbar", Nbars, sep = ".")), sep=":")

richness_PoissonBB_df <- tibble(estimate = unname(sapply(list_eb_fit_PoissonBB, function(x)
  total_richness(x)$lambda_post + ncol(x$feature_matrix)) ) ) %>%
  add_column(Model = "PoissonBB", 
             Nbar = rep(c("EB", Nbars), length(Ns)),
             n_train = rep(Ns, each = length(Nbars)+ 1)) 


richness_NegBinBB_df <- tibble(estimate = unname(sapply(list_eb_fit_NegBinBB, function(x)
  total_richness(x)$mu0_post + ncol(x$feature_matrix)) ) ) %>%
  add_column(Model = "NegBinBB", 
             Nbar = rep(c("EB", Nbars), length(Ns)),
             n_train = rep(Ns, each = length(Nbars)+ 1) ) 
# richness_NegBinBB_df <- tibble(estimate = c(unname(sapply(list_eb_fit_NegBinBB, function(x)
#   total_richness(x)$mu0_post + ncol(x$feature_matrix)) ), params_beta_NegBinBB$mu0) ) %>%
#   add_column(Model = "NegBinBB", 
#              Nbar = c(rep(c("EB", Nbars), length(Ns)), rep("E_emp",length(Ns))) ,
#              n_train = c(rep(Ns, each = length(Nbars)+ 1), Ns) ) 


joint_richness_long <- bind_rows(richness_PoissonBB_df, richness_NegBinBB_df) %>%
  mutate(Model = fct_relevel(Model, c( "NegBinBB", "PoissonBB")))

n_train.labs <- paste0("n = ", Ns,", Kn = ", Kn )
names(n_train.labs) <- Ns

# plots estimator
ggplot(joint_richness_long, aes( y=estimate, x=Model, shape = Nbar)) +
  geom_point(size = 2) +
  facet_wrap(.~ n_train, labeller = labeller(n_train = n_train.labs ), scales = "free_x", nrow = 1) +
  theme_light() +
  geom_hline(aes(yintercept = ncol(data_mat)), linetype = "dashed") +
  theme(legend.position = "top") +
  ylab("Posterior mean of N") +
  scale_y_continuous(breaks = pretty_breaks()) +
  rremove("xlab") +
  scale_shape_discrete(name = "Prior mean of N") +
  theme(aspect.ratio = 1)
ggsave(filename = "R_script_paper/Paper_plots/richness_points_beta_pis_eb.pdf", width = 8, height = 3.5, dpi = 300, units = "in", device='pdf')


# 2) Plot Richness: whole distributions for Nbar = 400 
list_eb_400_PoissonBB <- list_eb_fit_PoissonBB[grepl("Nbar.400", names(list_eb_fit_PoissonBB) )]
params_richness_400_PoissonBB <- tibble( lambda_prime = unname(sapply(list_eb_400_PoissonBB, function(x)
  total_richness(x)$lambda_post))) %>%
  add_column(n_train = Ns, Model = "PoissonBB") %>%
  mutate(lb = qpois(0.025, lambda_prime, lower.tail = TRUE, log.p = FALSE),
         ub = qpois(0.975, lambda_prime, lower.tail = TRUE, log.p = FALSE) )


list_eb_400_NegBinBB <- list_eb_fit_NegBinBB[grepl("Nbar.400", names(list_eb_fit_NegBinBB) )]
params_richness_400_NegBinBB <- tibble( n0_prime = unname(sapply(list_eb_400_NegBinBB, function(x)
  total_richness(x)$n0_post)),
  mu0_prime = unname(sapply(list_eb_400_NegBinBB, function(x)
    total_richness(x)$mu0_post))) %>%
  add_column(n_train = Ns, Model = "NegBinBB") %>%
  mutate(p_prime = 1/(mu0_prime/n0_prime + 1),
         lb = qnbinom(0.025, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE),
         ub = qnbinom(0.975, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE) )


bounds_distr <- tibble( lb = min(params_richness_400_PoissonBB$lb + Kn, 
                                 params_richness_400_NegBinBB$lb + Kn, H),
                        ub = max(params_richness_400_PoissonBB$ub + Kn, 
                                  params_richness_400_NegBinBB$ub + Kn, H))


dens_richness_PoissonBB_n15 <- tibble( x = bounds_distr$lb: bounds_distr$ub) %>%
  mutate( y = dpois(x - Kn[1] , lambda = params_richness_400_PoissonBB$lambda_prime[1])) %>%
  add_column(Model = "PoissonBB", n_train = params_richness_400_PoissonBB$n_train[1])

dens_richness_PoissonBB_n30 <- tibble( x = bounds_distr$lb : bounds_distr$ub) %>%
  mutate( y = dpois(x - Kn[2], lambda = params_richness_400_PoissonBB$lambda_prime[2])) %>%
  add_column(Model = "PoissonBB", n_train = params_richness_400_PoissonBB$n_train[2])

dens_richness_PoissonBB_n60 <- tibble( x = bounds_distr$lb : bounds_distr$ub) %>%
  mutate( y = dpois(x - Kn[3], lambda = params_richness_400_PoissonBB$lambda_prime[3])) %>%
  add_column(Model = "PoissonBB", n_train = params_richness_400_PoissonBB$n_train[3])


dens_richness_NegBinBB_n15 <- tibble( x = bounds_distr$lb : bounds_distr$ub) %>%
  mutate( y = dnbinom(x - Kn[1], size = params_richness_400_NegBinBB$n0_prime[1], 
                      prob = params_richness_400_NegBinBB$p_prime[1])) %>%
  add_column(Model = "NegBinBB", n_train = params_richness_400_NegBinBB$n_train[1])

dens_richness_NegBinBB_n30 <- tibble( x = bounds_distr$lb : bounds_distr$ub) %>%
  mutate( y = dnbinom(x - Kn[2], size = params_richness_400_NegBinBB$n0_prime[2], 
                      prob = params_richness_400_NegBinBB$p_prime[2])) %>%
  add_column(Model = "NegBinBB", n_train = params_richness_400_NegBinBB$n_train[2])

dens_richness_NegBinBB_n60 <- tibble( x = bounds_distr$lb : bounds_distr$ub) %>%
  mutate( y = dnbinom(x - Kn[3], size = params_richness_400_NegBinBB$n0_prime[3], 
                      prob = params_richness_400_NegBinBB$p_prime[3])) %>%
  add_column(Model = "NegBinBB", n_train = params_richness_400_NegBinBB$n_train[3]) 

dens_richnesses <- rbind(dens_richness_PoissonBB_n15,
                         dens_richness_PoissonBB_n30,
                         dens_richness_PoissonBB_n60,
                         dens_richness_NegBinBB_n15,
                         dens_richness_NegBinBB_n30,
                         dens_richness_NegBinBB_n60)




ggplot(dens_richnesses, aes(x = x, y = y, color = Model)) +
  geom_line() +
  geom_vline(aes(xintercept = sum(colSums(data_mat) > 0)), linetype="dashed") +
  facet_wrap(.~n_train, labeller = labeller(n_train = n_train.labs ), scales = "free_x") +
  theme_light() +
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  xlab("# distinct features") + rremove("ylab") +
  scale_color_tableau() +
  theme(aspect.ratio = 1)
ggsave(filename = "R_script_paper/Paper_plots/richness_Nbar400_beta_pis_eb.pdf", width = 8, height = 3.8, dpi = 300, units = "in", device='pdf')


# 3) Plot Extrapolation - EB version

# Extract accumulation curve of the observed sample (or average accumulation)
M = 300

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
list_eb_EB_PoissonBB <- list_eb_fit_PoissonBB[grepl("Nbar.emp", names(list_eb_fit_PoissonBB) )]
extr_PoissonBB_df <- tibble(lambda = unname(unlist(lapply(list_eb_EB_PoissonBB, function(x)
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

extr_PoissonBB_df$x <- as.integer(extr_PoissonBB_df$x)
extr_PoissonBB_df <- extr_PoissonBB_df %>%
  select(means, lb, ub, n_train, x, Model)

# NegBin
list_eb_EB_NegBinBB <- list_eb_fit_NegBinBB[grepl("Nbar.emp", names(list_eb_fit_NegBinBB) )]
extr_NegBinBB_df <- tibble(mu0 = unname(unlist(lapply(list_eb_EB_NegBinBB, function(x)
  extrapolation(object = x, M = M, seed = seed)$mu0_post ))),
  n0 = unname(unlist(lapply(list_eb_EB_NegBinBB, function(x)
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
             Model = "NegBinBB")
    
extr_NegBinBB_df$x <- as.integer(extr_NegBinBB_df$x)
extr_NegBinBB_df <- extr_NegBinBB_df %>%
  select(means, lb, ub, n_train, x, Model)


# GammaIBP
list_eb_EB_GammaIBP <- list_eb_fit_GammaIBP[grepl("Nbar.emp", names(list_eb_fit_GammaIBP) )]
extr_GammaIBP_df <- tibble(mu0 = unname(unlist(lapply(list_eb_EB_GammaIBP, function(x)
  extrapolation(object = x, M = M, seed = seed)$mu0_post ))),
  n0 = unname(unlist(lapply(list_eb_EB_GammaIBP, function(x)
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
             Model = "GammaIBP")

extr_GammaIBP_df$x <- as.integer(extr_GammaIBP_df$x)
extr_GammaIBP_df <- extr_GammaIBP_df %>%
  select(means, lb, ub, n_train, x, Model)

df_extr_GT_long <- list_extr_competitor_to_long(list_extr_GT, model = "GT") %>%
  rename(Model = model) %>%
  filter(t < n_train + M + 1)
df_extr_Chao_long <- list_extr_competitor_to_long(list_rare_extr_Chao, model = "Chao") %>%
  rename(Model = model) %>%
  filter(t < n_train + M + 1)

# Plot
temp <- tibble(n_train = Ns, xvalues = Ns)
extr_all_df <- rbind(extr_PoissonBB_df, extr_NegBinBB_df, extr_GammaIBP_df)

extr_all_df$Model <- factor(extr_all_df$Model, 
                            levels = c("PoissonBB", "NegBinBB", "GammaIBP",
                                       "Chao"))

ggplot(extr_all_df, aes(x, means, color = Model)) +
  geom_line(linetype = "dashed", linewidth = 0.9) +
  geom_point( data = accum_df_train, aes(x = x, y = n_feat),
              color="black", shape = 1, size = 1) +
  geom_point( data = accum_df_test, aes(x = x, y = n_feat),
              color="darkgrey", shape = 1, size = 0.1) +
  #geom_line(data = df_extr_GT_long, aes(t, value), linetype = "dashed", linewidth = 0.6) +
  geom_line(data = df_extr_Chao_long, aes(t, medians), linetype = "dashed", linewidth = 0.9) +
  facet_wrap(. ~ n_train,labeller = labeller(n_train = n_train.labs ),  scales = "free_x") +
  geom_vline(data = temp, mapping =  aes(xintercept = xvalues) , linetype = "dashed", color = "grey") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1)
ggsave(filename = "R_script_paper/Paper_plots/extr_beta_pis_eb.pdf", width = 8, height = 3.8, dpi = 300, units = "in", device='pdf')




### 1.C) MM approach ----

# Choose mechanism
mechanism = "beta_pis" # c("homogeneous", "random_uniform", "broken_stick", "log-normal")

# Fit and estimate richness, rarefaction and extrapolation for GibbsFA's (save workspace)
if (!file.exists(paste0("R_script_paper/eb_MM_",mechanism,"_fit_estimate_singledataset.RData"))) {
  
  c_fr <- 20
  
  # 3) GammaIBP
  
  # Initialization and known parameters
  eb_init_GammaIBP <- list(alpha = 0.5, s = 1, a = 1, b = 1)
  eb_known_GammaIBP <- list()
  
  # NegBin for different variances
  #cfrs <- c(2, 10, 50)
  
  # 4) Call the routine to perform simulations
  eb_MM_fit_estimate_ecological_scenario_singledataset(mechanism = mechanism, 
                                                       c_fr,
                                                       eb_init_GammaIBP, eb_known_GammaIBP, seed = 1234)
  
}

# Load the Work space
load(paste0("R_script_paper/eb_MM_",mechanism,"_fit_estimate_singledataset.RData"))

Kn <- sapply(Ns, function(n) sum(colSums(data_mat[1:n,]) > 0)  )
sum(colSums(data_mat) > 0)
plot(rarefaction(data_mat))

###### Comparison Poisson, NegBin and Gamma ------

# 0.B) Check on rarefaction
n_rare <- Ns[3]
eb_MM_biased_fit_PoissonBB_rare <- list_eb_MM_biased_fit_PoissonBB[[3]]
eb_MM_cens_fit_PoissonBB_rare <- list_eb_MM_cens_fit_PoissonBB[[3]]
eb_MM_biased_fit_NegBinBB_rare <- list_eb_MM_biased_fit_NegBinBB[[3]]
eb_MM_cens_fit_NegBinBB_rare <- list_eb_MM_cens_fit_NegBinBB[[3]]
eb_MM_fit_GammaIBP_rare <- list_eb_MM_fit_GammaIBP[[3]]

accum_df <- tibble( x = 0:n_rare,
                    n_feat = c(0,rarefaction(data_mat[1:n_rare,], n_reorderings = 3)))

rare_MM_biased_PoissonBB <- tibble( lambda_post = unname(unlist(
  rarefaction(object = eb_MM_biased_fit_PoissonBB_rare, seed = seed)$lambda_post ))) %>%
  rename(means = lambda_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "PoissonBB_MM_biased",
             x = c(1:n_rare,0))

rare_MM_cens_PoissonBB <- tibble( lambda_post = unname(unlist(
  rarefaction(object = eb_MM_cens_fit_PoissonBB_rare, seed = seed)$lambda_post ))) %>%
  rename(means = lambda_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "PoissonBB_MM_cens",
             x = c(1:n_rare,0))

rare_MM_biased_NegBinBB <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_MM_biased_fit_NegBinBB_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "NegBinBB_MM_biased",
             x = c(1:n_rare,0))

rare_MM_cens_NegBinBB <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_MM_cens_fit_NegBinBB_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "NegBinBB_MM_cens",
             x = c(1:n_rare,0))

rare_MM_GammaIBP <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_MM_fit_GammaIBP_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "GammaIBP_MM",
             x = c(1:n_rare,0)) 

df_rare <- rbind(#rare_MM_biased_PoissonBB, 
                 rare_MM_cens_PoissonBB,
                 #rare_MM_biased_NegBinBB, 
                 rare_MM_cens_NegBinBB
                 #rare_MM_GammaIBP)
                 )
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
ggsave(filename = paste0("R_script_paper/Paper_plots/rarefaction_",mechanism,"_eb_MM.pdf"), width = 8, height = 3.5, dpi = 300, units = "in", device='pdf')


# 0.c) Check on K_n_r
n_knr <- Ns[3]
eb_MM_biased_fit_PoissonBB_knr <- list_eb_MM_biased_fit_PoissonBB[[9]]
eb_MM_cens_fit_PoissonBB_knr <- list_eb_MM_cens_fit_PoissonBB[[9]]
eb_MM_biased_fit_NegBinBB_knr <- list_eb_MM_biased_fit_NegBinBB[[9]]
eb_MM_cens_fit_NegBinBB_knr <- list_eb_MM_cens_fit_NegBinBB[[9]]
eb_MM_fit_GammaIBP_knr <- list_eb_MM_fit_GammaIBP[[3]]


observed_K_n_r <- tibble( r = 1:n_knr,
                          k_n_r = K_n_r(data_mat[1:n_knr,], n_reorderings = 1)[[paste0('N = ', n_knr)]])

K_n_r_MM_biased_PoissonBB <- tibble( lambda_est = unname(unlist(
  K_n_r(object = eb_MM_biased_fit_PoissonBB_knr, n = n_knr)[[paste0('N = ', n_knr)]]$lambda_est ))) %>%
  rename(means = lambda_est) %>%
  add_column(Model = "PoissonBB_MM_biased",
             r = 1:n_knr)

K_n_r_MM_cens_PoissonBB <- tibble( lambda_est = unname(unlist(
  K_n_r(object = eb_MM_cens_fit_PoissonBB_knr, n = n_knr)[[paste0('N = ', n_knr)]]$lambda_est ))) %>%
  rename(means = lambda_est) %>%
  add_column(Model = "PoissonBB_MM_cens",
             r = 1:n_knr)

K_n_r_MM_biased_NegBinBB <- tibble( mu0_est = unname(unlist(
  K_n_r(object = eb_MM_biased_fit_NegBinBB_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
  rename(means = mu0_est) %>%
  add_column(Model = "NegBinBB_MM_biased",
             r = 1:n_knr)

K_n_r_MM_cens_NegBinBB <- tibble( mu0_est = unname(unlist(
  K_n_r(object = eb_MM_cens_fit_NegBinBB_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
  rename(means = mu0_est) %>%
  add_column(Model = "NegBinBB_MM_cens",
             r = 1:n_knr)

K_n_r_MM_GammaIBP <- tibble( mu0_est = unname(unlist(
  K_n_r(object = eb_MM_fit_GammaIBP_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
  rename(means = mu0_est) %>%
  add_column(Model = "GammaIBP_MM",
             r = 1:n_knr)

df_K_n_r <- rbind(K_n_r_MM_biased_PoissonBB, K_n_r_MM_cens_PoissonBB,
                  #K_n_r_MM_biased_NegBinBB, K_n_r_MM_cens_NegBinBB,
                  K_n_r_MM_GammaIBP)

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
ggsave(filename = paste0("R_script_paper/Paper_plots/knr_", mechanism, "_eb_MM.pdf"), width = 5, height = 4.5, dpi = 300, units = "in", device='pdf')



##### Prediction -----------

# 1) Plot Richness: Point plot expected value (number of features)
labels_comb_bb <- paste(rep(paste("n_train", Ns, sep = "."), each = length(Nbars)+1),
                        c("Nbar.emp" , paste("Nbar", Nbars, sep = ".")), sep=":")

richness_MM_biased_PoissonBB_df <- tibble(estimate = unname(sapply(list_eb_MM_biased_fit_PoissonBB, function(x)
  total_richness(x)$lambda_post + ncol(x$feature_matrix)) ) ) %>%
  add_column(Model = "PoissonBB_biased", 
             Nbar = rep(c("EB", Nbars), length(Ns)),
             n_train = rep(Ns, each = length(Nbars)+ 1)) 

richness_MM_cens_PoissonBB_df <- tibble(estimate = unname(sapply(list_eb_MM_cens_fit_PoissonBB, function(x)
  total_richness(x)$lambda_post + ncol(x$feature_matrix)) ) ) %>%
  add_column(Model = "PoissonBB_cens", 
             Nbar = rep(c("EB", Nbars), length(Ns)),
             n_train = rep(Ns, each = length(Nbars)+ 1)) 


richness_MM_biased_NegBinBB_df <- tibble(estimate = unname(sapply(list_eb_MM_biased_fit_NegBinBB, function(x)
  total_richness(x)$mu0_post + ncol(x$feature_matrix)) ) ) %>%
  add_column(Model = "NegBinBB_biased", 
             Nbar = rep(c("EB", Nbars), length(Ns)),
             n_train = rep(Ns, each = length(Nbars)+ 1) )

richness_MM_cens_NegBinBB_df <- tibble(estimate = unname(sapply(list_eb_MM_cens_fit_NegBinBB, function(x)
  total_richness(x)$mu0_post + ncol(x$feature_matrix)) ) ) %>%
  add_column(Model = "NegBinBB_cens", 
             Nbar = rep(c("EB", Nbars), length(Ns)),
             n_train = rep(Ns, each = length(Nbars)+ 1) )


joint_richness_long <- bind_rows(#richness_MM_biased_PoissonBB_df, 
                                 richness_MM_cens_PoissonBB_df,
                                 #richness_MM_biased_NegBinBB_df, 
                                 richness_MM_cens_NegBinBB_df) %>%
  mutate(Model = fct_relevel(Model))

n_train.labs <- paste0("n = ", Ns,", Kn = ", Kn )
names(n_train.labs) <- Ns

# plots estimator
ggplot(joint_richness_long, aes( y=estimate, x=Model, shape = Nbar)) +
  geom_point(size = 2) +
  facet_wrap(.~ n_train, labeller = labeller(n_train = n_train.labs ), scales = "free_x", nrow = 1) +
  theme_light() +
  geom_hline(aes(yintercept = ncol(data_mat)), linetype = "dashed") +
  theme(legend.position = "top") +
  ylab("Posterior mean of N") +
  scale_y_continuous(breaks = pretty_breaks()) +
  rremove("xlab") +
  scale_shape_discrete(name = "Prior mean of N") +
  theme(aspect.ratio = 1)
ggsave(filename = paste0("R_script_paper/Paper_plots/richness_points_", mechanism,"_eb_MM.pdf"), width = 8, height = 3.5, dpi = 300, units = "in", device='pdf')


# 2) Plot Richness: whole distributions for Nbar = 600 
list_eb_MM_biased_400_PoissonBB <- list_eb_MM_biased_fit_PoissonBB[grepl("Nbar.600", names(list_eb_MM_biased_fit_PoissonBB) )]
params_richness_MM_biased_400_PoissonBB <- tibble( lambda_prime = unname(sapply(list_eb_MM_biased_400_PoissonBB, function(x)
  total_richness(x)$lambda_post))) %>%
  add_column(n_train = Ns, Model = "PoissonBB_MM_biased") %>%
  mutate(lb = qpois(0.025, lambda_prime, lower.tail = TRUE, log.p = FALSE),
         ub = qpois(0.975, lambda_prime, lower.tail = TRUE, log.p = FALSE) )

list_eb_MM_cens_400_PoissonBB <- list_eb_MM_cens_fit_PoissonBB[grepl("Nbar.600", names(list_eb_MM_cens_fit_PoissonBB) )]
params_richness_MM_cens_400_PoissonBB <- tibble( lambda_prime = unname(sapply(list_eb_MM_cens_400_PoissonBB, function(x)
  total_richness(x)$lambda_post))) %>%
  add_column(n_train = Ns, Model = "PoissonBB_MM_cens") %>%
  mutate(lb = qpois(0.025, lambda_prime, lower.tail = TRUE, log.p = FALSE),
         ub = qpois(0.975, lambda_prime, lower.tail = TRUE, log.p = FALSE) )


list_eb_MM_biased_400_NegBinBB <- list_eb_MM_biased_fit_NegBinBB[grepl("Nbar.600", names(list_eb_MM_biased_fit_NegBinBB) )]
params_richness_MM_biased_400_NegBinBB <- tibble( n0_prime = unname(sapply(list_eb_MM_biased_400_NegBinBB, function(x)
  total_richness(x)$n0_post)),
  mu0_prime = unname(sapply(list_eb_MM_biased_400_NegBinBB, function(x)
    total_richness(x)$mu0_post))) %>%
  add_column(n_train = Ns, Model = "NegBinBB_MM_biased") %>%
  mutate(p_prime = 1/(mu0_prime/n0_prime + 1),
         lb = qnbinom(0.025, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE),
         ub = qnbinom(0.975, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE) )


list_eb_MM_cens_400_NegBinBB <- list_eb_MM_cens_fit_NegBinBB[grepl("Nbar.600", names(list_eb_MM_cens_fit_NegBinBB) )]
params_richness_MM_cens_400_NegBinBB <- tibble( n0_prime = unname(sapply(list_eb_MM_cens_400_NegBinBB, function(x)
  total_richness(x)$n0_post)),
  mu0_prime = unname(sapply(list_eb_MM_cens_400_NegBinBB, function(x)
    total_richness(x)$mu0_post))) %>%
  add_column(n_train = Ns, Model = "NegBinBB_MM_cens") %>%
  mutate(p_prime = 1/(mu0_prime/n0_prime + 1),
         lb = qnbinom(0.025, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE),
         ub = qnbinom(0.975, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE) )


bounds_distr <- tibble( lb = min(#params_richness_MM_biased_400_PoissonBB$lb + Kn,
                                 params_richness_MM_cens_400_PoissonBB$lb + Kn,
                                 #params_richness_MM_biased_400_NegBinBB$lb + Kn, 
                                 params_richness_MM_cens_400_NegBinBB$lb + Kn, 
                                 H),
                        ub = max(#params_richness_MM_biased_400_PoissonBB$ub + Kn,
                                 params_richness_MM_cens_400_PoissonBB$ub + Kn,
                                 #params_richness_MM_biased_400_NegBinBB$ub + Kn, 
                                 params_richness_MM_cens_400_NegBinBB$ub + Kn, 
                                 H))


dens_richness_MM_biased_PoissonBB <- bind_rows(lapply(1:length(Ns), function(j) 
  tibble( x = bounds_distr$lb: bounds_distr$ub) %>%
    mutate( y = dpois(x - Kn[j] , lambda = params_richness_MM_biased_400_PoissonBB$lambda_prime[j])) %>%
    add_column(Model = "PoissonBB_biased", n_train = params_richness_MM_biased_400_PoissonBB$n_train[j])
  ))

dens_richness_MM_cens_PoissonBB <- bind_rows(lapply(1:length(Ns), function(j) 
  tibble( x = bounds_distr$lb: bounds_distr$ub) %>%
    mutate( y = dpois(x - Kn[j] , lambda = params_richness_MM_cens_400_PoissonBB$lambda_prime[j])) %>%
    add_column(Model = "PoissonBB_cens", n_train = params_richness_MM_cens_400_PoissonBB$n_train[j])
))


# dens_richness_PoissonBB_n15 <- tibble( x = bounds_distr$lb: bounds_distr$ub) %>%
#   mutate( y = dpois(x - Kn[1] , lambda = params_richness_400_PoissonBB$lambda_prime[1])) %>%
#   add_column(Model = "PoissonBB", n_train = params_richness_400_PoissonBB$n_train[1])
# 
# dens_richness_PoissonBB_n30 <- tibble( x = bounds_distr$lb : bounds_distr$ub) %>%
#   mutate( y = dpois(x - Kn[2], lambda = params_richness_400_PoissonBB$lambda_prime[2])) %>%
#   add_column(Model = "PoissonBB", n_train = params_richness_400_PoissonBB$n_train[2])
# 
# dens_richness_PoissonBB_n60 <- tibble( x = bounds_distr$lb : bounds_distr$ub) %>%
#   mutate( y = dpois(x - Kn[3], lambda = params_richness_400_PoissonBB$lambda_prime[3])) %>%
#   add_column(Model = "PoissonBB", n_train = params_richness_400_PoissonBB$n_train[3])


dens_richness_MM_biased_NegBinBB <- bind_rows(lapply(1:length(Ns), function(j) 
  tibble( x = bounds_distr$lb: bounds_distr$ub) %>%
    mutate( y = dnbinom(x - Kn[j], size = params_richness_MM_biased_400_NegBinBB$n0_prime[j], 
                        prob = params_richness_MM_biased_400_NegBinBB$p_prime[j])) %>%
    add_column(Model = "NegBinBB_biased", n_train = params_richness_MM_biased_400_NegBinBB$n_train[j])
))

dens_richness_MM_cens_NegBinBB <- bind_rows(lapply(1:length(Ns), function(j) 
  tibble( x = bounds_distr$lb: bounds_distr$ub) %>%
    mutate( y = dnbinom(x - Kn[j], size = params_richness_MM_cens_400_NegBinBB$n0_prime[j], 
                        prob = params_richness_MM_cens_400_NegBinBB$p_prime[j])) %>%
    add_column(Model = "NegBinBB_cens", n_train = params_richness_MM_cens_400_NegBinBB$n_train[j])
))

# dens_richness_NegBinBB_n15 <- tibble( x = bounds_distr$lb : bounds_distr$ub) %>%
#   mutate( y = dnbinom(x - Kn[1], size = params_richness_400_NegBinBB$n0_prime[1], 
#                       prob = params_richness_400_NegBinBB$p_prime[1])) %>%
#   add_column(Model = "NegBinBB", n_train = params_richness_400_NegBinBB$n_train[1])
# 
# dens_richness_NegBinBB_n30 <- tibble( x = bounds_distr$lb : bounds_distr$ub) %>%
#   mutate( y = dnbinom(x - Kn[2], size = params_richness_400_NegBinBB$n0_prime[2], 
#                       prob = params_richness_400_NegBinBB$p_prime[2])) %>%
#   add_column(Model = "NegBinBB", n_train = params_richness_400_NegBinBB$n_train[2])
# 
# dens_richness_NegBinBB_n60 <- tibble( x = bounds_distr$lb : bounds_distr$ub) %>%
#   mutate( y = dnbinom(x - Kn[3], size = params_richness_400_NegBinBB$n0_prime[3], 
#                       prob = params_richness_400_NegBinBB$p_prime[3])) %>%
#   add_column(Model = "NegBinBB", n_train = params_richness_400_NegBinBB$n_train[3]) 

dens_richnesses <- rbind(#dens_richness_MM_biased_PoissonBB,
                         dens_richness_MM_cens_PoissonBB,
                         #dens_richness_MM_biased_NegBinBB,
                         dens_richness_MM_cens_NegBinBB)



n_train.labs <- paste0("n = ", Ns,", Kn = ", Kn )
names(n_train.labs) <- Ns

ggplot(dens_richnesses, aes(x = x, y = y, color = Model)) +
  geom_line() +
  geom_vline(aes(xintercept = H), linetype="dashed") +
  facet_wrap(.~n_train, labeller = labeller(n_train = n_train.labs ), scales = "free_x") +
  theme_light() +
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  xlab("# distinct features") + rremove("ylab") +
  scale_color_tableau() +
  theme(aspect.ratio = 1)
ggsave(filename = paste0("R_script_paper/Paper_plots/richness_Nbar600_", mechanism, "_eb_MM.pdf"), width = 8, height = 3.8, dpi = 300, units = "in", device='pdf')

lambda_MM_cens <- list_eb_MM_cens_400_PoissonBB$`n_train.120:Nbar.600`$lambda
alpha_MM_cens <- list_eb_MM_cens_400_PoissonBB$`n_train.120:Nbar.600`$alpha
theta_MM_cens <- list_eb_MM_cens_400_PoissonBB$`n_train.120:Nbar.600`$theta

var_N_prime_Poisson <- lambda_MM_cens*(1 - feature_fraction(120, alpha_MM_cens, theta_MM_cens))

mu0_MM_cens <- list_eb_MM_cens_400_NegBinBB$`n_train.120:Nbar.600`$mu0
n0_MM_cens <- list_eb_MM_cens_400_NegBinBB$`n_train.120:Nbar.600`$n0
p0_MM_cens <- 1/(mu0_MM_cens/n0_MM_cens + 1)

var_N_prime_NegBin <- mu0_MM_cens/p0_MM_cens



# 3) Plot Extrapolation - EB version

# Extract accumulation curve of the observed sample (or average accumulation)
M = 250

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
list_eb_MM_biased_PoissonBB <- list_eb_MM_biased_fit_PoissonBB[grepl("Nbar.emp", names(list_eb_MM_biased_fit_PoissonBB) )]
extr_MM_biased_PoissonBB_df <- tibble(lambda = unname(unlist(lapply(list_eb_MM_biased_PoissonBB, function(x)
  extrapolation(object = x, M = M, seed = seed)$lambda_post ))),
  n_train = rep(Ns, each = M),
  Kn = rep(Kn, each = M)) %>%
  mutate(lb = qpois(0.025, lambda, lower.tail = TRUE, log.p = FALSE),
         ub = qpois(0.975, lambda, lower.tail = TRUE, log.p = FALSE)) %>%
  rename(means = lambda) %>%
  add_row(means = 0, lb = 0, ub = 0, n_train = Ns, Kn = Kn) %>%
  mutate(means = means + Kn, lb = lb + Kn, ub = ub + Kn) %>%
  add_column(x = c(unlist(sapply(Ns, function(n) (n+1):(M+n))), Ns),
             Model = "PoissonBB_MM_biased")

extr_MM_biased_PoissonBB_df$x <- as.integer(extr_MM_biased_PoissonBB_df$x)
extr_MM_biased_PoissonBB_df <- extr_MM_biased_PoissonBB_df %>%
  select(means, lb, ub, n_train, x, Model)


list_eb_MM_cens_PoissonBB <- list_eb_MM_cens_fit_PoissonBB[grepl("Nbar.emp", names(list_eb_MM_cens_fit_PoissonBB) )]
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
list_eb_MM_biased_NegBinBB <- list_eb_MM_biased_fit_NegBinBB[grepl("Nbar.emp", names(list_eb_MM_biased_fit_NegBinBB) )]
extr_MM_biased_NegBinBB_df <- tibble(mu0 = unname(unlist(lapply(list_eb_MM_biased_NegBinBB, function(x)
  extrapolation(object = x, M = M, seed = seed)$mu0_post ))),
  n0 = unname(unlist(lapply(list_eb_MM_biased_NegBinBB, function(x)
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
             Model = "NegBinBB_MM_biased")

extr_MM_biased_NegBinBB_df$x <- as.integer(extr_MM_biased_NegBinBB_df$x)
extr_MM_biased_NegBinBB_df <- extr_MM_biased_NegBinBB_df %>%
  select(means, lb, ub, n_train, x, Model)


list_eb_MM_cens_NegBinBB <- list_eb_MM_cens_fit_NegBinBB[grepl("Nbar.emp", names(list_eb_MM_cens_fit_NegBinBB) )]
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
list_eb_MM_GammaIBP <- list_eb_MM_fit_GammaIBP[grepl("Nbar.emp", names(list_eb_MM_fit_GammaIBP) )]
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


# df_extr_GT_long <- list_extr_competitor_to_long(list_extr_GT, model = "GT") %>%
#   rename(Model = model) %>%
#   filter(t < n_train + M + 1)
# df_extr_Chao_long <- list_extr_competitor_to_long(list_rare_extr_Chao, model = "Chao") %>%
#   rename(Model = model) %>%
#   filter(t < n_train + M + 1)

# Plot
temp <- tibble(n_train = Ns, xvalues = Ns)
extr_all_df <- rbind(#extr_MM_biased_PoissonBB_df, 
  extr_MM_cens_PoissonBB_df,
  #extr_MM_biased_NegBinBB_df, 
  extr_MM_cens_NegBinBB_df,
  extr_MM_GammaIBP_df)

extr_all_df$Model <- factor(extr_all_df$Model)

ggplot(extr_all_df, aes(x, means, color = Model)) +
  geom_line(linetype = "dashed", linewidth = 0.9) +
  geom_ribbon(aes(ymin = lb, ymax = ub, color = Model), linewidth = 0.8, alpha = 0.1) +
  geom_point( data = accum_df_train, aes(x = x, y = n_feat),
              color="black", shape = 1, size = 1) +
  geom_point( data = accum_df_test, aes(x = x, y = n_feat),
              color="darkgrey", shape = 1, size = 0.1) +
  #geom_line(data = df_extr_GT_long, aes(t, value), linetype = "dashed", linewidth = 0.6) +
  #geom_line(data = df_extr_Chao_long, aes(t, medians), linetype = "dashed", linewidth = 0.9) +
  facet_wrap(. ~ n_train,labeller = labeller(n_train = n_train.labs ),  scales = "free_x") +
  geom_vline(data = temp, mapping =  aes(xintercept = xvalues) , linetype = "dashed", color = "grey") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1)
ggsave(filename = "R_script_paper/Paper_plots/extr_beta_pis_eb.pdf", width = 8, height = 3.8, dpi = 300, units = "in", device='pdf')



########### 2) Main script repeated-dataset -------


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


