
rm(list=ls())
library(ProductFormFA)
library(tidyverse)
library(scales)
library(ggpubr)
library(ggthemes)

source("R_script_paper/Routine_Chao.R")
source("R_script_paper/utils.R")


# Function to generate binary matrix according to 4 different mechanisms ----

generate_data <- function(mechanism, n, H, seed = 1234){

  if (!mechanism %in% c("homogeneous", "random_uniform", "broken_stick", "log-normal")){
    stop("Invalid generating mechanism.")
  }
  
  set.seed(seed)
  
  if (mechanism == "homogeneous"){ 
    
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
    
    # Fill the richness structures
    list_richness_PoissonBB[[lab_comb_bb]] <- total_richness(object = PoissonBB_fit)
    list_richness_NegBinBB[[lab_comb_bb]] <- total_richness(object = NegBinBB_fit)
    
    # Fill the rarefaction structures
    list_rare_PoissonBB[[lab_comb_bb]] <- rarefaction(object = PoissonBB_fit)
    list_rare_NegBinBB[[lab_comb_bb]] <- rarefaction(object = NegBinBB_fit)
    list_rare_GammaIBP[[lab_comb_ibp]] <- rarefaction(object = GammaIBP_fit)
    
    # Fill the extrapolation structures
    list_extr_PoissonBB[[lab_comb_bb]] <- extrapolation(object = PoissonBB_fit, M = M)
    list_extr_NegBinBB[[lab_comb_bb]] <- extrapolation(object = NegBinBB_fit, M = M)
    list_extr_GammaIBP[[lab_comb_ibp]] <- extrapolation(object = GammaIBP_fit, M = M)
    
    
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
      prior_obj_NegBinBB$mu0 <- 1/c_fr
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

# Choose mechanism
mechanism = "broken_stick" # c("homogeneous", "random_uniform", "broken_stick", "log-normal")
  
# Fit and estimate richness, rarefaction and extrapolation for GibbsFA's (save workspace)
if (!file.exists(paste0("R_script_paper/",mechanism,"_fit_estimate_singledataset.RData"))) {
  
  # 1) PoissonBB 
  
  # Initialization and MCMC setting 
  init_PoissonBB <- list(alpha_0 = -1, s_0 = 1)
  init_obj_PoissonBB <- initialization(model = "PoissonBB", init = init_PoissonBB )
  mcmcparams_PoissonBB <- list(tau = 0.1, S = 10000, n_burnin = 1000, thin = 2)
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
  add_column(Model= "Poisson",
             n_train = rep(Ns, each = length(list_richness_PoissonBB[[1]])*(length(Nbars)+1)),
             Nbar = rep(rep(c("emp", Nbars), each = length(list_richness_PoissonBB[[1]])), length(Ns) ) )

richness_NegBinBB_long <- gather(as_tibble(list_richness_NegBinBB), training, estimate,
                                labels_comb_bb,
                                factor_key=TRUE) %>%
  add_column(Model= "NegBin",
             n_train = rep(Ns, each = length(list_richness_PoissonBB[[1]])*(length(Nbars)+1)),
             Nbar = rep(rep(c("emp", Nbars), each = length(list_richness_PoissonBB[[1]])), length(Ns) ))

joint_richness_long <- bind_rows(richness_PoissonBB_long, richness_NegBinBB_long) %>%
  mutate(Model = fct_relevel(Model, c( "NegBin", "Poisson")))

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


# 2) Plot Richness: whole distributions for Empirical 
joint_emp_long <- joint_richness_long %>%
  filter(Nbar == "emp") %>%
  select(-training)

dev.new()
ggplot(joint_emp_long, aes(x = estimate, color = Model)) +
  stat_density(aes(x=estimate, colour=Model),
               geom="line",position="identity", adjust = 2.5) +
  geom_vline(aes(xintercept = 500), linetype="dashed") +
  facet_wrap(.~n_train, scales = "free_x") +
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


