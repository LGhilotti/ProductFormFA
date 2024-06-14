
rm(list=ls())
library(ProductFormFA)
library(tidyverse)
library(scales)
library(ggpubr)
library(ggthemes)

source("R_script_paper/Routine_Chao.R")
source("R_script_paper/utils.R")


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

eb_EFPF_fit_estimate_polynomial_scenario_singledataset <- function(xi, 
                                                                   vars_fct_NegBinBB,
                                                                   vars_GammaIBP,
                                                                   eb_init_BB, eb_known_BB,
                                                                   eb_init_IBP, eb_known_IBP,
                                                                   seed = 1234){
  
  if (xi <= 0){
    stop("Invalid generating mechanism.")
  }
  
  # Set training dimensions and dimension of the whole sample 
  Ns <- c(10, 50, 250) 
  N_max <- Ns[length(Ns)]
  L <- N_max + 400
  
  # Set maximum number of features
  H <- 10^6
  
  # Generate data 
  data_mat <- generate_data(xi = xi, n = L, H = H, seed = seed)
  data_mat <- data_mat[, colSums(data_mat) > 0]
  
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
    
  }
  
  # Save the entire workspace related to the mechanism just performed
  save(list = ls(all.names = TRUE), file =  paste0("R_script_paper/eb_EFPF_poly_",xi,"_fit_estimate_singledataset.RData"))
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



# 1) Main script single-dataset -------



##  EFPF approach ----

# Choose mechanism
xi = 1 

# Fit and estimate richness, rarefaction and extrapolation for GibbsFA's (save workspace)
if (!file.exists(paste0("R_script_paper/eb_EFPF_poly_",xi,"_fit_estimate_singledataset.RData"))) {
  
  
  vars_fct_NegBinBB <- c(10, 1000) 
  vars_GammaIBP <- c(0.01, 100)
  
  eb_init_BB <- list(alpha = -10, s = 10, Nhat_prime = 100)
  eb_known_BB <- list()
  
  eb_init_IBP <- list(alpha = 0.5, s = 1, Gamma = 10)
  eb_known_IBP <- list()
  
  # 4) Call the routine to perform simulations
  eb_EFPF_fit_estimate_polynomial_scenario_singledataset(xi = xi, 
                                                         vars_fct_NegBinBB,
                                                         vars_GammaIBP,
                                                         eb_init_BB, eb_known_BB,
                                                         eb_init_IBP, eb_known_IBP, 
                                                         seed = 123456)
  
}

# Load the Work space
load(paste0("R_script_paper/eb_EFPF_poly_",xi,"_fit_estimate_singledataset.RData"))

Kn <- sapply(Ns, function(n) sum(colSums(data_mat[1:n,]) > 0)  )
sum(colSums(data_mat) > 0)
plot(rarefaction(data_mat))

### Comparison between PoissonBB, NegBinBB, GammaIBP ----

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
                        levels = c("PoissonBB/NegBinBB", "GammaIBP"))

ggplot(accum_df, aes(x = x, y = n_feat)) +
  geom_point(color="black", shape = 18, size = 0.4) + 
  #facet_wrap(.~ Model, scales = "free_x", nrow = 1) +
  #geom_ribbon(aes(ymin = lb_bands, ymax = ub_bands), color = "red" , linewidth = 0.8, alpha = 0.1) +
  geom_line(data = df_rare, aes(x = x, y = means, color = Model), linetype = "dashed") + 
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) + 
  scale_color_tableau(
    labels = c(
      "PoissonBB/NegBinBB" = "Mixtures of BBs",
      "GammaIBP" = "Mixtures of IBPs"
    )
  ) 
ggsave(filename = paste0("R_script_paper/Paper_plots/rarefaction_poly_", xi, "_eb_EFPF.pdf"), width = 4, height = 4, dpi = 300, units = "in", device='pdf')


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
  filter(r < 10)

df_K_n_r_plot <- df_K_n_r %>%
  filter(r %in% c(r_positive$r))

observed_K_n_r_plot <- observed_K_n_r %>%
  filter(r %in% c(r_positive$r))

ggplot(observed_K_n_r_plot,  aes(x = r, y = k_n_r)) +
  geom_point(color="black", shape = 19, size = 0.5) +
  geom_line( data = df_K_n_r_plot, aes(x = r, y = means, color = Model), linetype = "dashed") +
  scale_y_log10() +
  xlab("r") + ylab(expression(m[r])) + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) + 
  scale_color_tableau(
    labels = c(
      "PoissonBB/NegBinBB" = "Mixtures of BBs",
      "GammaIBP" = "Mixtures of IBPs"
    )
  )
ggsave(filename = paste0("R_script_paper/Paper_plots/knr_poly_", xi, "_eb_EFPF.pdf"), width = 4, height = 4, dpi = 300, units = "in", device='pdf')




### Extrapolation ----

# GT estimators

if (!file.exists(paste0("R_script_paper/eb_Freq_poly_",xi,"_fit_estimate_singledataset.RData"))) {
  
  M <- 400
  
  # List to store smoothed GT's estimates
  list_extr_GT <- vector(mode="list")
  
  for (j in 1:length(Ns)){
    
    n_train <- Ns[j]
    
    train_mat <- data_mat[1:n_train,]
    train_mat <- train_mat[, colSums(train_mat) > 0]
    K <- ncol(train_mat)
    
    lab_comb_bb <- paste0("n_train.",n_train,":Nbar.emp")
    lab_comb_ibp <- paste("n_train", n_train, sep = ".")
    
    # Smoothed Good-Toulmin extrapolation
    
    # Compute SFS vector and CTS vector
    sfs <- tabulate(colSums(train_mat))
    cts <- sapply(2:n_train, function(n) ncol(train_mat[1:n,colSums(train_mat[1:n,]) > 0])   )
    cts <- c(0, sum(train_mat[1,]) , cts)
    
    list_extr_GT[[lab_comb_ibp]] <- predict_good_toulmin(n_train, M, sfs, cts, alternative = 0)$preds
    
    
  }
  
  # Save the entire workspace related to the mechanism just performed
  save(list = ls(all.names = TRUE), file =  paste0("R_script_paper/eb_Freq_poly_",xi,"_fit_estimate_singledataset.RData"))
  
}

# Load the Work space
load(paste0("R_script_paper/eb_Freq_poly_",xi,"_fit_estimate_singledataset.RData"))


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

# GammaIBP
extr_EFPF_GammaIBP_df <- tibble(means = numeric(), 
                                lb = numeric(), ub = numeric(),
                                n_train = integer(), 
                                x = integer(), Model = character())

for (var_GammaIBP in vars_GammaIBP){
  
  list_eb_EFPF_GammaIBP_var <- 
    list_eb_EFPF_fit_GammaIBP[[paste0("var.", var_GammaIBP)]]
  
  extr_EFPF_GammaIBP_df_var <- tibble(mu0 = unname(unlist(lapply(list_eb_EFPF_GammaIBP_var, function(x)
    extrapolation(object = x, M = M, seed = seed)$mu0_post ))),
    n0 = unname(unlist(lapply(list_eb_EFPF_GammaIBP_var, function(x)
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
               Model = paste0("GammaIBP, var = ", var_GammaIBP))
  
  extr_EFPF_GammaIBP_df_var$x <- as.integer(extr_EFPF_GammaIBP_df_var$x)
  extr_EFPF_GammaIBP_df_var <- extr_EFPF_GammaIBP_df_var %>%
    select(means, lb, ub, n_train, x, Model)
  
  
  extr_EFPF_GammaIBP_df <- bind_rows(extr_EFPF_GammaIBP_df, 
                                     extr_EFPF_GammaIBP_df_var)
  
}

df_extr_GT_long <- list_extr_competitor_to_long(list_extr_GT, model = "GT") %>%
  rename(Model = model) %>%
  filter(t < n_train + M + 1)

# Plot
temp <- tibble(n_train = Ns, xvalues = Ns)
extr_EFPF_PoissonBB_df <- extr_EFPF_PoissonBB_df %>%
  add_column(Model_gen = "PoissonBB/NegBinBB")
extr_EFPF_NegBinBB_df <- extr_EFPF_NegBinBB_df %>%
  add_column(Model_gen = "PoissonBB/NegBinBB")
extr_EFPF_GammaIBP_df <- extr_EFPF_GammaIBP_df %>%
  add_column(Model_gen = "GammaIBP")
extr_all_df <- rbind(#extr_EFPF_PoissonBB_df, 
                     #extr_EFPF_NegBinBB_df,
                     extr_EFPF_GammaIBP_df)

extr_all_df$Model <- factor(extr_all_df$Model, 
                            levels = c("PoissonBB", 
                                       paste0("NegBinBB x", vars_fct_NegBinBB),
                                       paste0("GammaIBP, var = ", vars_GammaIBP)))

n_train.labs <- paste0("n = ", Ns,", Kn = ", Kn )
names(n_train.labs) <- Ns

ggplot(extr_all_df, aes(x, means, color = Model)) +
  geom_line(linetype = "dashed", linewidth = 1) +
  geom_ribbon(aes(ymin = lb, ymax = ub, color = Model), linewidth = 0.8, alpha = 0) +
  geom_point( data = accum_df_train, aes(x = x, y = n_feat),
              color="black", shape = 21, size = 1) +
  geom_point( data = accum_df_test, aes(x = x, y = n_feat),
              color="black", shape = 21, size = 0.5) +
  geom_line(data = df_extr_GT_long, aes(t, value), linetype = "dashed", linewidth = 0.9) +
  facet_grid(. ~ n_train,  
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
ggsave(filename = paste0("R_script_paper/Paper_plots/extr_poly_", xi, "_eb_EFPF.pdf"), width = 9, height = 4, dpi = 300, units = "in", device='pdf')



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


