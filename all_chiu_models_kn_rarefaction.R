rm(list=ls())
library(ggmcmc)
library(coda)
library(ggpubr)
library(scales)
library(tidyverse)
library(gridExtra)
library(grid)
library(patchwork)
library(ProductFormFA)
library(ggthemes)

####
##### MODEL 1: the homogeneous model ############################
####

##### Single dataset -> Ntilde (Poiss/NB) and extrapolation (Poiss/NB/Gamma) #####

###### 1) Read results:  MCMC convergence ####################
load(file = "chiu_model_simulation/m3/m3_params_poiss.Rda")
load(file =  "chiu_model_simulation/m3/m3_params_negbin.Rda")
#load(file =  "chiu_model_simulation/m3/m3_params_negbin_prior.Rda")
load(file =  "chiu_model_simulation/m3/m3_params_ibp.Rda" )
#load(file =  "chiu_model_simulation/m3/m3_params_sp.Rda" )

###### 2) Read results: samples from limiting distributions (Poiss/NB) ##############
load(file = "chiu_model_simulation/m3/m3_ntilde_poiss.Rda")
load(file = "chiu_model_simulation/m3/m3_ntilde_negbin.Rda")
#load(file = "chiu_model_simulation/m3/m3_ntilde_negbin_prior.Rda")

###### 3) Read results: CI for extrapolation (Poiss/NB/Gamma) ################
list_kmn_pred_test_poiss <- readRDS(file = "chiu_model_simulation/m3/m3_ci_poiss.rds")
list_kmn_pred_test_negbin <- readRDS(file = "chiu_model_simulation/m3/m3_ci_negbin.rds")
#list_kmn_pred_test_negbin_prior <- readRDS(file = "chiu_model_simulation/m3/m3_ci_negbin_prior.rds")
list_kmn_pred_test_ibp <- readRDS(file = "chiu_model_simulation/m3/m3_ci_ibp.rds")
#list_kmn_pred_test_sp <- readRDS(file = "chiu_model_simulation/m3/m3_ci_sp.rds")

###### 4) Read the data ###############################
data_mat <- readRDS(file = "chiu_model_simulation/m3/m3_data_mat.rds")
data_list <- create_features_list(data_mat)
L <- nrow(data_mat)
Ms <- unique(sapply(list_kmn_pred_test_poiss, function(l) length(l$medians)))
Ns <- L - Ms

Nbars <- c(200,400,600)
c_fr <- 10


##### 9) Rarefaction curve and comparison with Chao bands #######

hors <- rep(L, length(Ns))

# Store objects
labels_comb_bb <- paste(rep(paste("N", Ns, sep = "."), each = length(Nbars)+1),
                        c(paste("Nbar", Nbars, sep = "."),"Nbar.emp"), sep=":")

labels_comb_ibp <- paste("N", Ns, sep = ".")

list_kn_rarefaction_poiss <- vector(mode="list", length = length(Ns)*(length(Nbars)+1))
names(list_kn_rarefaction_poiss) <- labels_comb_bb
list_kn_rarefaction_negbin <- vector(mode="list", length = length(Ns)*(length(Nbars)+1))
names(list_kn_rarefaction_negbin) <- labels_comb_bb
list_kn_rarefaction_ibp <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_ibp) <- labels_comb_ibp
# list_kn_rarefaction_sp <- vector(mode="list", length = length(Ns))
# names(list_kn_rarefaction_sp) <- paste("N", Ns, sep = ".")

for (j in 1:length(Ns)){
  N <- Ns[j]
  hor <- hors[j]
  M <- L - N
  
  train_mat <- data_mat[1:N,]
  test_mat <- data_mat[(N+1):L, ]
  # convert the binary matrix into list of features
  train_list <- create_features_list(train_mat)
  test_list <- create_features_list(test_mat)
  
  ##### 2.2) Run the empirical Nbar 
  Nbar_emp <- beta_binomial_estimator(train_mat)
  
  # Set prior hyperparameters specific for poisson, with fixed lambda
  lambda_poiss <- Nbar_emp
  
  # Set prior hyperparameters specific for NB, with fixed parameters
  nstar_nb <- Nbar_emp/(c_fr - 1)
  p_nb <- 1/c_fr
  
  # Label for accessing element of structures related to N and Nbar
  lab_comb_bb <- paste0("N.",N,":Nbar.emp")
  lab_comb_ibp <- paste0("N.",N)
  
  lab_alpha_bb <- paste0("alpha:",lab_comb_bb)
  lab_theta_bb <- paste0("theta:",lab_comb_bb)
  lab_alpha_ibp <- paste0("alpha:",lab_comb_ibp)
  lab_theta_ibp <- paste0("theta:",lab_comb_ibp)
  lab_a <- paste0("a:",lab_comb_ibp)
  lab_b <- paste0("b:",lab_comb_ibp)
  
  # Poisson
  alpha_chain_poiss <- params_poiss[[lab_alpha_bb]] 
  theta_chain_poiss <- params_poiss[[lab_theta_bb]] 
  
  kn_chain_poiss <- generate_Kmn_chain_poiss(lambda_poiss, alpha_chain_poiss,
                                             theta_chain_poiss, M = hor, n=0)
  
  est_ci_poiss <- matrix(NA, nrow = hor, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:hor){
    est_ci_poiss[m,] <- quantile(kn_chain_poiss[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_poiss <- list("medians" = est_ci_poiss[,2],
                       "lbs" = est_ci_poiss[,1],
                       "ubs" = est_ci_poiss[,3])
  
  list_kn_rarefaction_poiss[[lab_comb_bb]] <- est_ci_poiss
  
  # NegBin
  alpha_chain_negbin <- params_negbin[[lab_alpha_bb]] 
  theta_chain_negbin <- params_negbin[[lab_theta_bb]]
  
  kn_chain_negbin <- generate_Kmn_chain_negbin(nstar_nb, p_nb,
                                               alpha_chain_negbin, theta_chain_negbin,
                                               M = hor, n=0)
  
  est_ci_negbin <- matrix(NA, nrow = hor, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:hor){
    est_ci_negbin[m,] <- quantile(kn_chain_negbin[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_negbin <- list("medians" = est_ci_negbin[,2],
                        "lbs" = est_ci_negbin[,1],
                        "ubs" = est_ci_negbin[,3])
  
  list_kn_rarefaction_negbin[[lab_comb_bb]] <- est_ci_negbin
  
  # IBP
  a_chain_ibp <- params_ibp[[lab_a]]
  b_chain_ibp <- params_ibp[[lab_b]] 
  alpha_chain_ibp <- params_ibp[[lab_alpha_ibp]] 
  theta_chain_ibp <- params_ibp[[lab_theta_ibp]] 
  
  kn_chain_ibp <- generate_Kmn_chain_gamma_ibp(a_chain_ibp, b_chain_ibp, alpha_chain_ibp, 
                                               theta_chain_ibp, M = hor, n=0)
  
  est_ci_ibp <- matrix(NA, nrow = hor, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:hor){
    est_ci_ibp[m,] <- quantile(kn_chain_ibp[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_ibp <- list("medians" = est_ci_ibp[,2],
                     "lbs" = est_ci_ibp[,1],
                     "ubs" = est_ci_ibp[,3])
  
  list_kn_rarefaction_ibp[[lab_comb_ibp]] <- est_ci_ibp
  
  # # SB-SP
  # c_chain_sp <- params_sp[[paste0("c.",N)]]
  # beta_chain_sp <- params_sp[[paste0("beta.",N)]] 
  # alpha_chain_sp <- params_sp[[paste0("alpha.",N)]] 
  # 
  # kn_chain_sp <- generate_Kmn_chain_gamma_ibp(a_chain = c_chain_sp + 1,
  #                                             b_chain = beta_chain_sp*(1-alpha_chain_sp)/alpha_chain_sp, 
  #                                             alpha_chain = alpha_chain_sp,
  #                                             theta_chain = 1 - alpha_chain_sp,
  #                                             M = N, n = 0)
  # 
  # est_ci_sp <- matrix(NA, nrow = N, ncol = 3)
  # # first column = lower bound
  # # second columns = medians
  # # third columns = upper bound
  # for (m in 1:N){
  #   est_ci_sp[m,] <- quantile(kn_chain_sp[m,], probs = c(0.025,0.5,0.975))
  # }
  # est_ci_sp <- list("medians" = est_ci_sp[,2],
  #                    "lbs" = est_ci_sp[,1],
  #                    "ubs" = est_ci_sp[,3])
  # 
  # list_kn_rarefaction_sp[[paste0("N.",N)]] <- est_ci_sp
  # 
  
  #### 2.3) Run the Nbars on the grid ###
  for (v in 1:length(Nbars)){
    
    Nbar <- Nbars[v]
    
    # Set prior hyperparameters specific for poisson, with fixed lambda
    lambda_poiss <- Nbar
    
    # Set prior hyperparameters specific for NB, with fixed parameters
    nstar_nb <- Nbar/(c_fr - 1)
    p_nb <- 1/c_fr
    
    # Label for accessing element of structures related to N and Nbar
    lab_comb <- paste0("N.",N,":Nbar.",Nbar)
    
    lab_alpha <- paste0("alpha:",lab_comb)
    lab_theta <- paste0("theta:",lab_comb)
    
    # Poisson
    alpha_chain_poiss <- params_poiss[[lab_alpha]] 
    theta_chain_poiss <- params_poiss[[lab_theta]] 
    
    kn_chain_poiss <- generate_Kmn_chain_poiss(lambda_poiss, alpha_chain_poiss,
                                               theta_chain_poiss, M = hor, n=0)
    
    est_ci_poiss <- matrix(NA, nrow = hor, ncol = 3)
    # first column = lower bound
    # second columns = medians
    # third columns = upper bound
    for (m in 1:hor){
      est_ci_poiss[m,] <- quantile(kn_chain_poiss[m,], probs = c(0.025,0.5,0.975))
    }
    est_ci_poiss <- list("medians" = est_ci_poiss[,2],
                         "lbs" = est_ci_poiss[,1],
                         "ubs" = est_ci_poiss[,3])
    
    list_kn_rarefaction_poiss[[lab_comb]] <- est_ci_poiss
    
    # NegBin
    alpha_chain_negbin <- params_negbin[[lab_alpha]] 
    theta_chain_negbin <- params_negbin[[lab_theta]]
    
    kn_chain_negbin <- generate_Kmn_chain_negbin(nstar_nb, p_nb,
                                                 alpha_chain_negbin, theta_chain_negbin,
                                                 M = hor, n=0)
    
    est_ci_negbin <- matrix(NA, nrow = hor, ncol = 3)
    # first column = lower bound
    # second columns = medians
    # third columns = upper bound
    for (m in 1:hor){
      est_ci_negbin[m,] <- quantile(kn_chain_negbin[m,], probs = c(0.025,0.5,0.975))
    }
    est_ci_negbin <- list("medians" = est_ci_negbin[,2],
                          "lbs" = est_ci_negbin[,1],
                          "ubs" = est_ci_negbin[,3])
    
    list_kn_rarefaction_negbin[[lab_comb]] <- est_ci_negbin
    
  }
}

# # Poisson
# saveRDS(list_kn_rarefaction_poiss, "chiu_model_simulation/m3/m3_ci_rare_poiss.rds")
# # Negative Binomial
# saveRDS(list_kn_rarefaction_negbin, "chiu_model_simulation/m3/m3_ci_rare_negbin.rds")
# # Gamma IBP
# saveRDS(list_kn_rarefaction_ibp, "chiu_model_simulation/m3/m3_ci_rare_ibp.rds")


#### 9.b) Chao's bands for rarefaction/extrapolation ####
source("Chao_code/Subroutine_for_iNEXT.R")

list_chao_rare <- vector(mode="list", length = length(Ns))
names(list_chao_rare) <- labels_comb_ibp

for (j in 1:length(Ns)){
  
  N <- Ns[j]
  hor <- hors[j]
  
  train_mat <- data_mat[1:N,]
  # convert the binary matrix into list of features
  train_list <- create_features_list(train_mat)
  
  lab_comb <- paste0("N.",N)
  
  # Determine the frequency vector of the training sets
  Q_vec <- colSums(train_mat)
  Q_vec <- Q_vec[Q_vec>0]
  
  # Compute the curves with confidence intervals
  chao_res <- iNEXT.Sam(Spec = Q_vec, T = N, endpoint = hor)
  
  chao_res_rare <- as_tibble(chao_res[["q=0"]]) %>%
    select(-Cov.hat) %>%
    rename(medians = D0.hat, lbs = Norm.CI.Low, ubs = Norm.CI.High)
  
  
  list_chao_rare[[lab_comb]] <- as.data.frame(chao_res_rare)

}
# Richness
saveRDS(list_chao_rare, "chiu_model_simulation/m3/m3_chao_rare.rds")

