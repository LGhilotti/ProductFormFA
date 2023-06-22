####
##### MODEL 1: the homogeneous model ############################
####

rm(list=ls())
library(ggmcmc)
library(coda)
library(ggpubr)
library(scales)
library(tidyverse)

##### Single dataset -> Ntilde (Poiss/NB) and extrapolation (Poiss/NB/Gamma) #####

###### 1) Read results ####################
load(file = "chao_model_simulation/m1/m1_params_poiss.Rda")
load(file =  "chao_model_simulation/m1/m1_params_negbin.Rda")
load(file =  "chao_model_simulation/m1/m1_params_ibp.Rda" )
list_kmn_pred_test_poiss <- readRDS(file = "chao_model_simulation/m1/m1_ci_poiss.rds")

###### 4) Read the data ###############################
data_mat <- readRDS(file = "chao_model_simulation/m1/m1_data_mat.rds")
L <- nrow(data_mat)
Ms <- sapply(list_kmn_pred_test_poiss, function(l) length(l$medians))
Ns <- L - Ms

list_kn_rarefaction_poiss <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_poiss) <- paste("N", Ns, sep = ".")
list_kn_rarefaction_negbin <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_negbin) <- paste("N", Ns, sep = ".")
list_kn_rarefaction_ibp <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_ibp) <- paste("N", Ns, sep = ".")

for (j in 1:length(Ns)){
  N <- Ns[j]
  M <- L - N
  
  train_mat <- data_mat[1:N,]
  test_mat <- data_mat[(N+1):L, ]
  # convert the binary matrix into list of features
  train_list <- create_features_list(train_mat)
  test_list <- create_features_list(test_mat)
  
  # Poisson
  lambda_chain_poiss <- params_poiss[[paste0("lambda.",N)]]
  alpha_chain_poiss <- params_poiss[[paste0("alpha.",N)]] 
  theta_chain_poiss <- params_poiss[[paste0("theta.",N)]] 
  
  kn_chain_poiss <- generate_Kmn_chain_poiss(lambda_chain_poiss, alpha_chain_poiss,
                                             theta_chain_poiss, M = N, n=0)
  
  est_ci_poiss <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_poiss[m,] <- quantile(kn_chain_poiss[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_poiss <- list("medians" = est_ci_poiss[,2],
                       "lbs" = est_ci_poiss[,1],
                       "ubs" = est_ci_poiss[,3])
  
  list_kn_rarefaction_poiss[[paste0("N.",N)]] <- est_ci_poiss
  
  # NegBin
  nstar_chain_negbin <- params_negbin[[paste0("nstar.",N)]]
  p_chain_negbin <- params_negbin[[paste0("p.",N)]]
  alpha_chain_negbin <- params_negbin[[paste0("alpha.",N)]] 
  theta_chain_negbin <- params_negbin[[paste0("theta.",N)]]
  
  kn_chain_negbin <- generate_Kmn_chain_negbin(nstar_chain_negbin, p_chain_negbin,
                                               alpha_chain_negbin, theta_chain_negbin,
                                               M = N, n=0)
  
  est_ci_negbin <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_negbin[m,] <- quantile(kn_chain_negbin[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_negbin <- list("medians" = est_ci_negbin[,2],
                        "lbs" = est_ci_negbin[,1],
                        "ubs" = est_ci_negbin[,3])
  
  list_kn_rarefaction_negbin[[paste0("N.",N)]] <- est_ci_negbin
  
  # IBP
  a_chain_ibp <- params_ibp[[paste0("a.",N)]]
  b_chain_ibp <- params_ibp[[paste0("b.",N)]] 
  alpha_chain_ibp <- params_ibp[[paste0("alpha.",N)]] 
  theta_chain_ibp <- params_ibp[[paste0("theta.",N)]] 
  
  kn_chain_ibp <- generate_Kmn_chain_gamma_ibp(a_chain_ibp, b_chain_ibp, alpha_chain_ibp, 
                                               theta_chain_ibp, M = N, n=0)
  
  est_ci_ibp <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_ibp[m,] <- quantile(kn_chain_ibp[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_ibp <- list("medians" = est_ci_ibp[,2],
                     "lbs" = est_ci_ibp[,1],
                     "ubs" = est_ci_ibp[,3])
  
  list_kn_rarefaction_ibp[[paste0("N.",N)]] <- est_ci_ibp
  
 
}
  
 
# Poisson
saveRDS(list_kn_rarefaction_poiss, "chao_model_simulation/m1/m1_ci_insample_poiss.rds")
# Negative Binomial
saveRDS(list_kn_rarefaction_negbin, "chao_model_simulation/m1/m1_ci_insample_negbin.rds")
# Gamma IBP
saveRDS(list_kn_rarefaction_ibp, "chao_model_simulation/m1/m1_ci_insample_ibp.rds")

  

####
##### MODEL 2: the random uniform model ############################
####

rm(list=ls())
library(ggmcmc)
library(coda)
library(ggpubr)
library(scales)

##### Single dataset -> Ntilde (Poiss/NB) and extrapolation (Poiss/NB/Gamma) #####

###### 1) Read results ####################
load(file = "chao_model_simulation/m2/m2_params_poiss.Rda")
load(file =  "chao_model_simulation/m2/m2_params_negbin.Rda")
load(file =  "chao_model_simulation/m2/m2_params_ibp.Rda", )
list_kmn_pred_test_poiss <- readRDS(file = "chao_model_simulation/m2/m2_ci_poiss.rds")

###### 4) Read the data ###############################
data_mat <- readRDS(file = "chao_model_simulation/m2/m2_data_mat.rds")
L <- nrow(data_mat)
Ms <- sapply(list_kmn_pred_test_poiss, function(l) length(l$medians))
Ns <- L - Ms

list_kn_rarefaction_poiss <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_poiss) <- paste("N", Ns, sep = ".")
list_kn_rarefaction_negbin <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_negbin) <- paste("N", Ns, sep = ".")
list_kn_rarefaction_ibp <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_ibp) <- paste("N", Ns, sep = ".")

for (j in 1:length(Ns)){
  N <- Ns[j]
  M <- L - N
  
  train_mat <- data_mat[1:N,]
  test_mat <- data_mat[(N+1):L, ]
  # convert the binary matrix into list of features
  train_list <- create_features_list(train_mat)
  test_list <- create_features_list(test_mat)
  
  # Poisson
  lambda_chain_poiss <- params_poiss[[paste0("lambda.",N)]]
  alpha_chain_poiss <- params_poiss[[paste0("alpha.",N)]] 
  theta_chain_poiss <- params_poiss[[paste0("theta.",N)]] 
  
  kn_chain_poiss <- generate_Kmn_chain_poiss(lambda_chain_poiss, alpha_chain_poiss,
                                             theta_chain_poiss, M = N, n=0)
  
  est_ci_poiss <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_poiss[m,] <- quantile(kn_chain_poiss[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_poiss <- list("medians" = est_ci_poiss[,2],
                       "lbs" = est_ci_poiss[,1],
                       "ubs" = est_ci_poiss[,3])
  
  list_kn_rarefaction_poiss[[paste0("N.",N)]] <- est_ci_poiss
  
  # NegBin
  nstar_chain_negbin <- params_negbin[[paste0("nstar.",N)]]
  p_chain_negbin <- params_negbin[[paste0("p.",N)]]
  alpha_chain_negbin <- params_negbin[[paste0("alpha.",N)]] 
  theta_chain_negbin <- params_negbin[[paste0("theta.",N)]]
  
  kn_chain_negbin <- generate_Kmn_chain_negbin(nstar_chain_negbin, p_chain_negbin,
                                               alpha_chain_negbin, theta_chain_negbin,
                                               M = N, n=0)
  
  est_ci_negbin <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_negbin[m,] <- quantile(kn_chain_negbin[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_negbin <- list("medians" = est_ci_negbin[,2],
                        "lbs" = est_ci_negbin[,1],
                        "ubs" = est_ci_negbin[,3])
  
  list_kn_rarefaction_negbin[[paste0("N.",N)]] <- est_ci_negbin
  
  # IBP
  a_chain_ibp <- params_ibp[[paste0("a.",N)]]
  b_chain_ibp <- params_ibp[[paste0("b.",N)]] 
  alpha_chain_ibp <- params_ibp[[paste0("alpha.",N)]] 
  theta_chain_ibp <- params_ibp[[paste0("theta.",N)]] 
  
  kn_chain_ibp <- generate_Kmn_chain_gamma_ibp(a_chain_ibp, b_chain_ibp, alpha_chain_ibp, 
                                               theta_chain_ibp, M = N, n=0)
  
  est_ci_ibp <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_ibp[m,] <- quantile(kn_chain_ibp[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_ibp <- list("medians" = est_ci_ibp[,2],
                     "lbs" = est_ci_ibp[,1],
                     "ubs" = est_ci_ibp[,3])
  
  list_kn_rarefaction_ibp[[paste0("N.",N)]] <- est_ci_ibp
  
  
}


# Poisson
saveRDS(list_kn_rarefaction_poiss, "chao_model_simulation/m2/m2_ci_insample_poiss.rds")
# Negative Binomial
saveRDS(list_kn_rarefaction_negbin, "chao_model_simulation/m2/m2_ci_insample_negbin.rds")
# Gamma IBP
saveRDS(list_kn_rarefaction_ibp, "chao_model_simulation/m2/m2_ci_insample_ibp.rds")



####
##### MODEL 3: the broken stick model ############################
####

rm(list=ls())
library(ggmcmc)
library(coda)
library(ggpubr)
library(scales)

##### Single dataset -> Ntilde (Poiss/NB) and extrapolation (Poiss/NB/Gamma) #####

###### 1) Read results ####################
load(file = "chao_model_simulation/m3/m3_params_poiss.Rda")
load(file =  "chao_model_simulation/m3/m3_params_negbin.Rda")
load(file =  "chao_model_simulation/m3/m3_params_ibp.Rda", )
list_kmn_pred_test_poiss <- readRDS(file = "chao_model_simulation/m3/m3_ci_poiss.rds")

###### 4) Read the data ###############################
data_mat <- readRDS(file = "chao_model_simulation/m3/m3_data_mat.rds")
L <- nrow(data_mat)
Ms <- sapply(list_kmn_pred_test_poiss, function(l) length(l$medians))
Ns <- L - Ms

list_kn_rarefaction_poiss <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_poiss) <- paste("N", Ns, sep = ".")
list_kn_rarefaction_negbin <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_negbin) <- paste("N", Ns, sep = ".")
list_kn_rarefaction_ibp <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_ibp) <- paste("N", Ns, sep = ".")

for (j in 1:length(Ns)){
  N <- Ns[j]
  M <- L - N
  
  train_mat <- data_mat[1:N,]
  test_mat <- data_mat[(N+1):L, ]
  # convert the binary matrix into list of features
  train_list <- create_features_list(train_mat)
  test_list <- create_features_list(test_mat)
  
  # Poisson
  lambda_chain_poiss <- params_poiss[[paste0("lambda.",N)]]
  alpha_chain_poiss <- params_poiss[[paste0("alpha.",N)]] 
  theta_chain_poiss <- params_poiss[[paste0("theta.",N)]] 
  
  kn_chain_poiss <- generate_Kmn_chain_poiss(lambda_chain_poiss, alpha_chain_poiss,
                                             theta_chain_poiss, M = N, n=0)
  
  est_ci_poiss <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_poiss[m,] <- quantile(kn_chain_poiss[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_poiss <- list("medians" = est_ci_poiss[,2],
                       "lbs" = est_ci_poiss[,1],
                       "ubs" = est_ci_poiss[,3])
  
  list_kn_rarefaction_poiss[[paste0("N.",N)]] <- est_ci_poiss
  
  # NegBin
  nstar_chain_negbin <- params_negbin[[paste0("nstar.",N)]]
  p_chain_negbin <- params_negbin[[paste0("p.",N)]]
  alpha_chain_negbin <- params_negbin[[paste0("alpha.",N)]] 
  theta_chain_negbin <- params_negbin[[paste0("theta.",N)]]
  
  kn_chain_negbin <- generate_Kmn_chain_negbin(nstar_chain_negbin, p_chain_negbin,
                                               alpha_chain_negbin, theta_chain_negbin,
                                               M = N, n=0)
  
  est_ci_negbin <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_negbin[m,] <- quantile(kn_chain_negbin[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_negbin <- list("medians" = est_ci_negbin[,2],
                        "lbs" = est_ci_negbin[,1],
                        "ubs" = est_ci_negbin[,3])
  
  list_kn_rarefaction_negbin[[paste0("N.",N)]] <- est_ci_negbin
  
  # IBP
  a_chain_ibp <- params_ibp[[paste0("a.",N)]]
  b_chain_ibp <- params_ibp[[paste0("b.",N)]] 
  alpha_chain_ibp <- params_ibp[[paste0("alpha.",N)]] 
  theta_chain_ibp <- params_ibp[[paste0("theta.",N)]] 
  
  kn_chain_ibp <- generate_Kmn_chain_gamma_ibp(a_chain_ibp, b_chain_ibp, alpha_chain_ibp, 
                                               theta_chain_ibp, M = N, n=0)
  
  est_ci_ibp <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_ibp[m,] <- quantile(kn_chain_ibp[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_ibp <- list("medians" = est_ci_ibp[,2],
                     "lbs" = est_ci_ibp[,1],
                     "ubs" = est_ci_ibp[,3])
  
  list_kn_rarefaction_ibp[[paste0("N.",N)]] <- est_ci_ibp
  
  
}


# Poisson
saveRDS(list_kn_rarefaction_poiss, "chao_model_simulation/m3/m3_ci_insample_poiss.rds")
# Negative Binomial
saveRDS(list_kn_rarefaction_negbin, "chao_model_simulation/m3/m3_ci_insample_negbin.rds")
# Gamma IBP
saveRDS(list_kn_rarefaction_ibp, "chao_model_simulation/m3/m3_ci_insample_ibp.rds")


####
##### MODEL 4: the log-normal model ############################
####

rm(list=ls())
library(ggmcmc)
library(coda)
library(ggpubr)
library(scales)

##### Single dataset -> Ntilde (Poiss/NB) and extrapolation (Poiss/NB/Gamma) #####

###### 1) Read results ####################
load(file = "chao_model_simulation/m4/m4_params_poiss.Rda")
load(file =  "chao_model_simulation/m4/m4_params_negbin.Rda")
load(file =  "chao_model_simulation/m4/m4_params_ibp.Rda", )
list_kmn_pred_test_poiss <- readRDS(file = "chao_model_simulation/m4/m4_ci_poiss.rds")

###### 4) Read the data ###############################
data_mat <- readRDS(file = "chao_model_simulation/m4/m4_data_mat.rds")
L <- nrow(data_mat)
Ms <- sapply(list_kmn_pred_test_poiss, function(l) length(l$medians))
Ns <- L - Ms

list_kn_rarefaction_poiss <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_poiss) <- paste("N", Ns, sep = ".")
list_kn_rarefaction_negbin <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_negbin) <- paste("N", Ns, sep = ".")
list_kn_rarefaction_ibp <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_ibp) <- paste("N", Ns, sep = ".")

for (j in 1:length(Ns)){
  N <- Ns[j]
  M <- L - N
  
  train_mat <- data_mat[1:N,]
  test_mat <- data_mat[(N+1):L, ]
  # convert the binary matrix into list of features
  train_list <- create_features_list(train_mat)
  test_list <- create_features_list(test_mat)
  
  # Poisson
  lambda_chain_poiss <- params_poiss[[paste0("lambda.",N)]]
  alpha_chain_poiss <- params_poiss[[paste0("alpha.",N)]] 
  theta_chain_poiss <- params_poiss[[paste0("theta.",N)]] 
  
  kn_chain_poiss <- generate_Kmn_chain_poiss(lambda_chain_poiss, alpha_chain_poiss,
                                             theta_chain_poiss, M = N, n=0)
  
  est_ci_poiss <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_poiss[m,] <- quantile(kn_chain_poiss[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_poiss <- list("medians" = est_ci_poiss[,2],
                       "lbs" = est_ci_poiss[,1],
                       "ubs" = est_ci_poiss[,3])
  
  list_kn_rarefaction_poiss[[paste0("N.",N)]] <- est_ci_poiss
  
  # NegBin
  nstar_chain_negbin <- params_negbin[[paste0("nstar.",N)]]
  p_chain_negbin <- params_negbin[[paste0("p.",N)]]
  alpha_chain_negbin <- params_negbin[[paste0("alpha.",N)]] 
  theta_chain_negbin <- params_negbin[[paste0("theta.",N)]]
  
  kn_chain_negbin <- generate_Kmn_chain_negbin(nstar_chain_negbin, p_chain_negbin,
                                               alpha_chain_negbin, theta_chain_negbin,
                                               M = N, n=0)
  
  est_ci_negbin <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_negbin[m,] <- quantile(kn_chain_negbin[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_negbin <- list("medians" = est_ci_negbin[,2],
                        "lbs" = est_ci_negbin[,1],
                        "ubs" = est_ci_negbin[,3])
  
  list_kn_rarefaction_negbin[[paste0("N.",N)]] <- est_ci_negbin
  
  # IBP
  a_chain_ibp <- params_ibp[[paste0("a.",N)]]
  b_chain_ibp <- params_ibp[[paste0("b.",N)]] 
  alpha_chain_ibp <- params_ibp[[paste0("alpha.",N)]] 
  theta_chain_ibp <- params_ibp[[paste0("theta.",N)]] 
  
  kn_chain_ibp <- generate_Kmn_chain_gamma_ibp(a_chain_ibp, b_chain_ibp, alpha_chain_ibp, 
                                               theta_chain_ibp, M = N, n=0)
  
  est_ci_ibp <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_ibp[m,] <- quantile(kn_chain_ibp[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_ibp <- list("medians" = est_ci_ibp[,2],
                     "lbs" = est_ci_ibp[,1],
                     "ubs" = est_ci_ibp[,3])
  
  list_kn_rarefaction_ibp[[paste0("N.",N)]] <- est_ci_ibp
  
  
}


# Poisson
saveRDS(list_kn_rarefaction_poiss, "chao_model_simulation/m4/m4_ci_insample_poiss.rds")
# Negative Binomial
saveRDS(list_kn_rarefaction_negbin, "chao_model_simulation/m4/m4_ci_insample_negbin.rds")
# Gamma IBP
saveRDS(list_kn_rarefaction_ibp, "chao_model_simulation/m4/m4_ci_insample_ibp.rds")



####
##### MODEL 5: the Zipf-Mandelbrot model ############################
####

rm(list=ls())
library(ggmcmc)
library(coda)
library(ggpubr)
library(scales)

##### Single dataset -> Ntilde (Poiss/NB) and extrapolation (Poiss/NB/Gamma) #####

###### 1) Read results ####################
load(file = "chao_model_simulation/m5/m5_params_poiss.Rda")
load(file =  "chao_model_simulation/m5/m5_params_negbin.Rda")
load(file =  "chao_model_simulation/m5/m5_params_ibp.Rda" )
list_kmn_pred_test_poiss <- readRDS(file = "chao_model_simulation/m5/m5_ci_poiss.rds")

###### 4) Read the data ###############################
data_mat <- readRDS(file = "chao_model_simulation/m5/m5_data_mat.rds")
L <- nrow(data_mat)
Ms <- sapply(list_kmn_pred_test_poiss, function(l) length(l$medians))
Ns <- L - Ms

list_kn_rarefaction_poiss <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_poiss) <- paste("N", Ns, sep = ".")
list_kn_rarefaction_negbin <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_negbin) <- paste("N", Ns, sep = ".")
list_kn_rarefaction_ibp <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_ibp) <- paste("N", Ns, sep = ".")

for (j in 1:length(Ns)){
  N <- Ns[j]
  M <- L - N
  
  train_mat <- data_mat[1:N,]
  test_mat <- data_mat[(N+1):L, ]
  # convert the binary matrix into list of features
  train_list <- create_features_list(train_mat)
  test_list <- create_features_list(test_mat)
  
  # Poisson
  lambda_chain_poiss <- params_poiss[[paste0("lambda.",N)]]
  alpha_chain_poiss <- params_poiss[[paste0("alpha.",N)]] 
  theta_chain_poiss <- params_poiss[[paste0("theta.",N)]] 
  
  kn_chain_poiss <- generate_Kmn_chain_poiss(lambda_chain_poiss, alpha_chain_poiss,
                                             theta_chain_poiss, M = N, n=0)
  
  est_ci_poiss <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_poiss[m,] <- quantile(kn_chain_poiss[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_poiss <- list("medians" = est_ci_poiss[,2],
                       "lbs" = est_ci_poiss[,1],
                       "ubs" = est_ci_poiss[,3])
  
  list_kn_rarefaction_poiss[[paste0("N.",N)]] <- est_ci_poiss
  
  # NegBin
  nstar_chain_negbin <- params_negbin[[paste0("nstar.",N)]]
  p_chain_negbin <- params_negbin[[paste0("p.",N)]]
  alpha_chain_negbin <- params_negbin[[paste0("alpha.",N)]] 
  theta_chain_negbin <- params_negbin[[paste0("theta.",N)]]
  
  kn_chain_negbin <- generate_Kmn_chain_negbin(nstar_chain_negbin, p_chain_negbin,
                                               alpha_chain_negbin, theta_chain_negbin,
                                               M = N, n=0)
  
  est_ci_negbin <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_negbin[m,] <- quantile(kn_chain_negbin[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_negbin <- list("medians" = est_ci_negbin[,2],
                        "lbs" = est_ci_negbin[,1],
                        "ubs" = est_ci_negbin[,3])
  
  list_kn_rarefaction_negbin[[paste0("N.",N)]] <- est_ci_negbin
  
  # IBP
  a_chain_ibp <- params_ibp[[paste0("a.",N)]]
  b_chain_ibp <- params_ibp[[paste0("b.",N)]] 
  alpha_chain_ibp <- params_ibp[[paste0("alpha.",N)]] 
  theta_chain_ibp <- params_ibp[[paste0("theta.",N)]] 
  
  kn_chain_ibp <- generate_Kmn_chain_gamma_ibp(a_chain_ibp, b_chain_ibp, alpha_chain_ibp, 
                                               theta_chain_ibp, M = N, n=0)
  
  est_ci_ibp <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_ibp[m,] <- quantile(kn_chain_ibp[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_ibp <- list("medians" = est_ci_ibp[,2],
                     "lbs" = est_ci_ibp[,1],
                     "ubs" = est_ci_ibp[,3])
  
  list_kn_rarefaction_ibp[[paste0("N.",N)]] <- est_ci_ibp
  
  
}


# Poisson
saveRDS(list_kn_rarefaction_poiss, "chao_model_simulation/m5/m5_ci_insample_poiss.rds")
# Negative Binomial
saveRDS(list_kn_rarefaction_negbin, "chao_model_simulation/m5/m5_ci_insample_negbin.rds")
# Gamma IBP
saveRDS(list_kn_rarefaction_ibp, "chao_model_simulation/m5/m5_ci_insample_ibp.rds")




####
##### Unbounded-features scenario: Polynomial (exponent: 1) ############################
####

rm(list=ls())
library(ggmcmc)
library(coda)
library(ggpubr)
library(scales)

##### Single dataset -> Ntilde (Poiss/NB) and extrapolation (Poiss/NB/Gamma) #####

###### 1) Read results ####################
load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_params_poiss.Rda")
load(file =  "unbounded_features_simulation/unb_poly_1/unb_poly_1_params_negbin.Rda")
load(file =  "unbounded_features_simulation/unb_poly_1/unb_poly_1_params_ibp.Rda", )
list_kmn_pred_test_poiss <- readRDS(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_ci_poiss.rds")

###### 4) Read the data ###############################
data_mat <- readRDS(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_data_mat.rds")
L <- nrow(data_mat)
Ms <- sapply(list_kmn_pred_test_poiss, function(l) length(l$medians))
Ns <- L - Ms

list_kn_rarefaction_poiss <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_poiss) <- paste("N", Ns, sep = ".")
list_kn_rarefaction_negbin <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_negbin) <- paste("N", Ns, sep = ".")
list_kn_rarefaction_ibp <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_ibp) <- paste("N", Ns, sep = ".")

for (j in 1:length(Ns)){
  N <- Ns[j]
  M <- L - N
  
  train_mat <- data_mat[1:N,]
  test_mat <- data_mat[(N+1):L, ]
  # convert the binary matrix into list of features
  train_list <- create_features_list(train_mat)
  test_list <- create_features_list(test_mat)
  
  # Poisson
  lambda_chain_poiss <- params_poiss[[paste0("lambda.",N)]]
  alpha_chain_poiss <- params_poiss[[paste0("alpha.",N)]] 
  theta_chain_poiss <- params_poiss[[paste0("theta.",N)]] 
  
  kn_chain_poiss <- generate_Kmn_chain_poiss(lambda_chain_poiss, alpha_chain_poiss,
                                             theta_chain_poiss, M = N, n=0)
  
  est_ci_poiss <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_poiss[m,] <- quantile(kn_chain_poiss[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_poiss <- list("medians" = est_ci_poiss[,2],
                       "lbs" = est_ci_poiss[,1],
                       "ubs" = est_ci_poiss[,3])
  
  list_kn_rarefaction_poiss[[paste0("N.",N)]] <- est_ci_poiss
  
  # NegBin
  nstar_chain_negbin <- params_negbin[[paste0("nstar.",N)]]
  p_chain_negbin <- params_negbin[[paste0("p.",N)]]
  alpha_chain_negbin <- params_negbin[[paste0("alpha.",N)]] 
  theta_chain_negbin <- params_negbin[[paste0("theta.",N)]]
  
  kn_chain_negbin <- generate_Kmn_chain_negbin(nstar_chain_negbin, p_chain_negbin,
                                               alpha_chain_negbin, theta_chain_negbin,
                                               M = N, n=0)
  
  est_ci_negbin <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_negbin[m,] <- quantile(kn_chain_negbin[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_negbin <- list("medians" = est_ci_negbin[,2],
                        "lbs" = est_ci_negbin[,1],
                        "ubs" = est_ci_negbin[,3])
  
  list_kn_rarefaction_negbin[[paste0("N.",N)]] <- est_ci_negbin
  
  # IBP
  a_chain_ibp <- params_ibp[[paste0("a.",N)]]
  b_chain_ibp <- params_ibp[[paste0("b.",N)]] 
  alpha_chain_ibp <- params_ibp[[paste0("alpha.",N)]] 
  theta_chain_ibp <- params_ibp[[paste0("theta.",N)]] 
  
  kn_chain_ibp <- generate_Kmn_chain_gamma_ibp(a_chain_ibp, b_chain_ibp, alpha_chain_ibp, 
                                               theta_chain_ibp, M = N, n=0)
  
  est_ci_ibp <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_ibp[m,] <- quantile(kn_chain_ibp[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_ibp <- list("medians" = est_ci_ibp[,2],
                     "lbs" = est_ci_ibp[,1],
                     "ubs" = est_ci_ibp[,3])
  
  list_kn_rarefaction_ibp[[paste0("N.",N)]] <- est_ci_ibp
  
  
}


# Poisson
saveRDS(list_kn_rarefaction_poiss, "unbounded_features_simulation/unb_poly_1/unb_poly_1_ci_insample_poiss.rds")
# Negative Binomial
saveRDS(list_kn_rarefaction_negbin, "unbounded_features_simulation/unb_poly_1/unb_poly_1_ci_insample_negbin.rds")
# Gamma IBP
saveRDS(list_kn_rarefaction_ibp, "unbounded_features_simulation/unb_poly_1/unb_poly_1_ci_insample_ibp.rds")


####
##### Unbounded-features scenario: Polynomial (exponent: 1.2) ############################
####

rm(list=ls())
library(ggmcmc)
library(coda)
library(ggpubr)
library(scales)

##### Single dataset -> Ntilde (Poiss/NB) and extrapolation (Poiss/NB/Gamma) #####

###### 1) Read results ####################
load(file = "unbounded_features_simulation/unb_poly_1_2/unb_poly_1_2_params_poiss.Rda")
load(file =  "unbounded_features_simulation/unb_poly_1_2/unb_poly_1_2_params_negbin.Rda")
load(file =  "unbounded_features_simulation/unb_poly_1_2/unb_poly_1_2_params_ibp.Rda", )
list_kmn_pred_test_poiss <- readRDS(file = "unbounded_features_simulation/unb_poly_1_2/unb_poly_1_2_ci_poiss.rds")

###### 4) Read the data ###############################
data_mat <- readRDS(file = "unbounded_features_simulation/unb_poly_1_2/unb_poly_1_2_data_mat.rds")
L <- nrow(data_mat)
Ms <- sapply(list_kmn_pred_test_poiss, function(l) length(l$medians))
Ns <- L - Ms

list_kn_rarefaction_poiss <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_poiss) <- paste("N", Ns, sep = ".")
list_kn_rarefaction_negbin <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_negbin) <- paste("N", Ns, sep = ".")
list_kn_rarefaction_ibp <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_ibp) <- paste("N", Ns, sep = ".")

for (j in 1:length(Ns)){
  N <- Ns[j]
  M <- L - N
  
  train_mat <- data_mat[1:N,]
  test_mat <- data_mat[(N+1):L, ]
  # convert the binary matrix into list of features
  train_list <- create_features_list(train_mat)
  test_list <- create_features_list(test_mat)
  
  # Poisson
  lambda_chain_poiss <- params_poiss[[paste0("lambda.",N)]]
  alpha_chain_poiss <- params_poiss[[paste0("alpha.",N)]] 
  theta_chain_poiss <- params_poiss[[paste0("theta.",N)]] 
  
  kn_chain_poiss <- generate_Kmn_chain_poiss(lambda_chain_poiss, alpha_chain_poiss,
                                             theta_chain_poiss, M = N, n=0)
  
  est_ci_poiss <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_poiss[m,] <- quantile(kn_chain_poiss[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_poiss <- list("medians" = est_ci_poiss[,2],
                       "lbs" = est_ci_poiss[,1],
                       "ubs" = est_ci_poiss[,3])
  
  list_kn_rarefaction_poiss[[paste0("N.",N)]] <- est_ci_poiss
  
  # NegBin
  nstar_chain_negbin <- params_negbin[[paste0("nstar.",N)]]
  p_chain_negbin <- params_negbin[[paste0("p.",N)]]
  alpha_chain_negbin <- params_negbin[[paste0("alpha.",N)]] 
  theta_chain_negbin <- params_negbin[[paste0("theta.",N)]]
  
  kn_chain_negbin <- generate_Kmn_chain_negbin(nstar_chain_negbin, p_chain_negbin,
                                               alpha_chain_negbin, theta_chain_negbin,
                                               M = N, n=0)
  
  est_ci_negbin <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_negbin[m,] <- quantile(kn_chain_negbin[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_negbin <- list("medians" = est_ci_negbin[,2],
                        "lbs" = est_ci_negbin[,1],
                        "ubs" = est_ci_negbin[,3])
  
  list_kn_rarefaction_negbin[[paste0("N.",N)]] <- est_ci_negbin
  
  # IBP
  a_chain_ibp <- params_ibp[[paste0("a.",N)]]
  b_chain_ibp <- params_ibp[[paste0("b.",N)]] 
  alpha_chain_ibp <- params_ibp[[paste0("alpha.",N)]] 
  theta_chain_ibp <- params_ibp[[paste0("theta.",N)]] 
  
  kn_chain_ibp <- generate_Kmn_chain_gamma_ibp(a_chain_ibp, b_chain_ibp, alpha_chain_ibp, 
                                               theta_chain_ibp, M = N, n=0)
  
  est_ci_ibp <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_ibp[m,] <- quantile(kn_chain_ibp[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_ibp <- list("medians" = est_ci_ibp[,2],
                     "lbs" = est_ci_ibp[,1],
                     "ubs" = est_ci_ibp[,3])
  
  list_kn_rarefaction_ibp[[paste0("N.",N)]] <- est_ci_ibp
  
  
}


# Poisson
saveRDS(list_kn_rarefaction_poiss, "unbounded_features_simulation/unb_poly_1_2/unb_poly_1_2_ci_insample_poiss.rds")
# Negative Binomial
saveRDS(list_kn_rarefaction_negbin, "unbounded_features_simulation/unb_poly_1_2/unb_poly_1_2_ci_insample_negbin.rds")
# Gamma IBP
saveRDS(list_kn_rarefaction_ibp, "unbounded_features_simulation/unb_poly_1_2/unb_poly_1_2_ci_insample_ibp.rds")


####
##### Unbounded-features scenario: Polynomial (exponent: 0.8) ############################
####

rm(list=ls())
library(ggmcmc)
library(coda)
library(ggpubr)
library(scales)

##### Single dataset -> Ntilde (Poiss/NB) and extrapolation (Poiss/NB/Gamma) #####

###### 1) Read results ####################
load(file = "unbounded_features_simulation/unb_poly_0_8/unb_poly_0_8_params_poiss.Rda")
load(file =  "unbounded_features_simulation/unb_poly_0_8/unb_poly_0_8_params_negbin.Rda")
load(file =  "unbounded_features_simulation/unb_poly_0_8/unb_poly_0_8_params_ibp.Rda", )
list_kmn_pred_test_poiss <- readRDS(file = "unbounded_features_simulation/unb_poly_0_8/unb_poly_0_8_ci_poiss.rds")

###### 4) Read the data ###############################
data_mat <- readRDS(file = "unbounded_features_simulation/unb_poly_0_8/unb_poly_0_8_data_mat.rds")
L <- nrow(data_mat)
Ms <- sapply(list_kmn_pred_test_poiss, function(l) length(l$medians))
Ns <- L - Ms

list_kn_rarefaction_poiss <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_poiss) <- paste("N", Ns, sep = ".")
list_kn_rarefaction_negbin <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_negbin) <- paste("N", Ns, sep = ".")
list_kn_rarefaction_ibp <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_ibp) <- paste("N", Ns, sep = ".")

for (j in 1:length(Ns)){
  N <- Ns[j]
  M <- L - N
  
  train_mat <- data_mat[1:N,]
  test_mat <- data_mat[(N+1):L, ]
  # convert the binary matrix into list of features
  train_list <- create_features_list(train_mat)
  test_list <- create_features_list(test_mat)
  
  # Poisson
  lambda_chain_poiss <- params_poiss[[paste0("lambda.",N)]]
  alpha_chain_poiss <- params_poiss[[paste0("alpha.",N)]] 
  theta_chain_poiss <- params_poiss[[paste0("theta.",N)]] 
  
  kn_chain_poiss <- generate_Kmn_chain_poiss(lambda_chain_poiss, alpha_chain_poiss,
                                             theta_chain_poiss, M = N, n=0)
  
  est_ci_poiss <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_poiss[m,] <- quantile(kn_chain_poiss[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_poiss <- list("medians" = est_ci_poiss[,2],
                       "lbs" = est_ci_poiss[,1],
                       "ubs" = est_ci_poiss[,3])
  
  list_kn_rarefaction_poiss[[paste0("N.",N)]] <- est_ci_poiss
  
  # NegBin
  nstar_chain_negbin <- params_negbin[[paste0("nstar.",N)]]
  p_chain_negbin <- params_negbin[[paste0("p.",N)]]
  alpha_chain_negbin <- params_negbin[[paste0("alpha.",N)]] 
  theta_chain_negbin <- params_negbin[[paste0("theta.",N)]]
  
  kn_chain_negbin <- generate_Kmn_chain_negbin(nstar_chain_negbin, p_chain_negbin,
                                               alpha_chain_negbin, theta_chain_negbin,
                                               M = N, n=0)
  
  est_ci_negbin <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_negbin[m,] <- quantile(kn_chain_negbin[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_negbin <- list("medians" = est_ci_negbin[,2],
                        "lbs" = est_ci_negbin[,1],
                        "ubs" = est_ci_negbin[,3])
  
  list_kn_rarefaction_negbin[[paste0("N.",N)]] <- est_ci_negbin
  
  # IBP
  a_chain_ibp <- params_ibp[[paste0("a.",N)]]
  b_chain_ibp <- params_ibp[[paste0("b.",N)]] 
  alpha_chain_ibp <- params_ibp[[paste0("alpha.",N)]] 
  theta_chain_ibp <- params_ibp[[paste0("theta.",N)]] 
  
  kn_chain_ibp <- generate_Kmn_chain_gamma_ibp(a_chain_ibp, b_chain_ibp, alpha_chain_ibp, 
                                               theta_chain_ibp, M = N, n=0)
  
  est_ci_ibp <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_ibp[m,] <- quantile(kn_chain_ibp[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_ibp <- list("medians" = est_ci_ibp[,2],
                     "lbs" = est_ci_ibp[,1],
                     "ubs" = est_ci_ibp[,3])
  
  list_kn_rarefaction_ibp[[paste0("N.",N)]] <- est_ci_ibp
  
  
}


# Poisson
saveRDS(list_kn_rarefaction_poiss, "unbounded_features_simulation/unb_poly_0_8/unb_poly_0_8_ci_insample_poiss.rds")
# Negative Binomial
saveRDS(list_kn_rarefaction_negbin, "unbounded_features_simulation/unb_poly_0_8/unb_poly_0_8_ci_insample_negbin.rds")
# Gamma IBP
saveRDS(list_kn_rarefaction_ibp, "unbounded_features_simulation/unb_poly_0_8/unb_poly_0_8_ci_insample_ibp.rds")
