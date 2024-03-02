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
##### Fungi data ############################
####

###### 1) Read results:  MCMC convergence ####################
load(file = "mazziotta2016_application/fungi/mazz_fungi_op_params_poiss.Rda")
load(file =  "mazziotta2016_application/fungi/mazz_fungi_op_params_negbin.Rda")
load(file =  "mazziotta2016_application/fungi/mazz_fungi_op_params_ibp.Rda" )

###### 3) Read results: CI for extrapolation (Poiss/NB/Gamma) ################
list_kmn_pred_test_poiss <- readRDS(file = "mazziotta2016_application/fungi/mazz_fungi_op_ci_pred_poiss.rds")
list_kmn_pred_test_negbin <- readRDS(file = "mazziotta2016_application/fungi/mazz_fungi_op_ci_pred_negbin.rds")
list_kmn_pred_test_ibp <- readRDS(file = "mazziotta2016_application/fungi/mazz_fungi_op_ci_pred_ibp.rds")

###### 4) Read the data ###############################
data_mat <- readRDS(file = "mazziotta2016_application/fungi/mazz_fungi_data_mat.rds")
data_list <- create_features_list(data_mat)
L <- nrow(data_mat)

# Set prior hyperparameters specific for poisson, with fixed lambda
Nbar_emp <- beta_binomial_estimator(data_mat)

# Set prior hyperparameters specific for poisson, with fixed lambda
lambda_poiss <- Nbar_emp

# Set prior hyperparameters specific for NB, with fixed parameters
c_fr <- 10
nstar_nb <- Nbar_emp/(c_fr - 1)
p_nb <- 1/c_fr

###### 2) Run the algorithms ###############
labels_comb <- "Nbar.emp"
number_saved_iterations_poiss <- length(params_poiss$`alpha:Nbar.emp`)
number_saved_iterations_negbin <- length(params_negbin$`alpha:Nbar.emp`)
number_saved_iterations_ibp <- length(params_ibp$`alpha:Nbar.emp`)

# Choose the horizon for extrapolation
hor <- 1000

list_kn_rare_poiss <- vector(mode="list", length = 1)
names(list_kn_rare_poiss) <- labels_comb
list_kn_rare_negbin <- vector(mode="list", length = 1)
names(list_kn_rare_negbin) <- labels_comb
list_kn_rare_ibp <- vector(mode="list", length = 1)
names(list_kn_rare_ibp) <- labels_comb

list_kmn_pred_poiss <- vector(mode="list", length = 1)
names(list_kmn_pred_poiss) <- labels_comb
list_kmn_pred_negbin <- vector(mode="list", length = 1)
names(list_kmn_pred_negbin) <- labels_comb
list_kmn_pred_ibp <- vector(mode="list", length = 1)
names(list_kmn_pred_ibp) <- labels_comb


Kn = ncol(data_mat[,colSums(data_mat) > 0])
print(paste0("Number of observed features: ", Kn))

# Label for accessing element of structures related to N and Nbar

lab_alpha <- paste0("alpha:",labels_comb)
lab_theta <- paste0("theta:",labels_comb)
lab_a <- paste0("a:",labels_comb)
lab_b <- paste0("b:",labels_comb)
lab_c <- paste0("c:",labels_comb)
lab_beta <- paste0("beta:",labels_comb)


##### 9) Rarefaction curve #######

# Poisson
alpha_chain_poiss <- params_poiss[[lab_alpha]] 
theta_chain_poiss <- params_poiss[[lab_theta]] 

kn_chain_poiss <- generate_Kmn_chain_poiss(lambda_poiss, alpha_chain_poiss,
                                           theta_chain_poiss, M = L , n=0)

est_ci_poiss <- matrix(NA, nrow = L, ncol = 3)
# first column = lower bound
# second columns = medians
# third columns = upper bound
for (m in 1:L){
  est_ci_poiss[m,c(1,3)] <- quantile(kn_chain_poiss[m,], probs = c(0.025,0.975))
  est_ci_poiss[m, 2] <- mean(kn_chain_poiss[m,])
}
est_ci_poiss <- list("means" = est_ci_poiss[,2],
                     "lbs" = est_ci_poiss[,1],
                     "ubs" = est_ci_poiss[,3])

list_kn_rare_poiss[[labels_comb]] <- est_ci_poiss

# NegBin
alpha_chain_negbin <- params_negbin[[lab_alpha]] 
theta_chain_negbin <- params_negbin[[lab_theta]]

kn_chain_negbin <- generate_Kmn_chain_negbin(nstar_nb, p_nb,
                                             alpha_chain_negbin, theta_chain_negbin,
                                             M = L, n=0)

est_ci_negbin <- matrix(NA, nrow = L, ncol = 3)
# first column = lower bound
# second columns = medians
# third columns = upper bound
for (m in 1:L){
  est_ci_negbin[m,c(1,3)] <- quantile(kn_chain_negbin[m,], probs = c(0.025,0.975))
  est_ci_negbin[m, 2] <- mean(kn_chain_negbin[m,])
}
est_ci_negbin <- list("means" = est_ci_negbin[,2],
                      "lbs" = est_ci_negbin[,1],
                      "ubs" = est_ci_negbin[,3])

list_kn_rare_negbin[[labels_comb]] <- est_ci_negbin

# # IBP
# a_chain_ibp <- params_ibp[[lab_a]]
# b_chain_ibp <- params_ibp[[lab_b]] 
# alpha_chain_ibp <- params_ibp[[lab_alpha]] 
# theta_chain_ibp <- params_ibp[[lab_theta]] 
# 
# kn_chain_ibp <- generate_Kmn_chain_gamma_ibp(a_chain_ibp, b_chain_ibp, alpha_chain_ibp, 
#                                              theta_chain_ibp, M = L, n=0)
# 
# est_ci_ibp <- matrix(NA, nrow = L, ncol = 3)
# # first column = lower bound
# # second columns = medians
# # third columns = upper bound
# for (m in 1:L){
#   est_ci_ibp[m,c(1,3)] <- quantile(kn_chain_ibp[m,], probs = c(0.025,0.975))
#   est_ci_ibp[m, 2] <- mean(kn_chain_ibp[m,])
# }
# est_ci_ibp <- list("means" = est_ci_ibp[,2],
#                    "lbs" = est_ci_ibp[,1],
#                    "ubs" = est_ci_ibp[,3])
# 
# list_kn_rare_ibp[[labels_comb]] <- est_ci_ibp




# Poisson
saveRDS(list_kn_rare_poiss, "mazziotta2016_application/fungi/mazz_fungi_op_ci_rare_poiss.rds")
# Negative Binomial
saveRDS(list_kn_rare_negbin, "mazziotta2016_application/fungi/mazz_fungi_op_ci_rare_negbin.rds")
# Gamma IBP
#saveRDS(list_kn_rare_ibp, "mazziotta2016_application/fungi/mazz_fungi_op_ci_rare_ibp.rds")


######## 7) Extrapolation in the future set (Poiss/NB/Gamma) ----------

# Poisson
kmn_chain_poiss <- generate_Kmn_chain_poiss(lambda_poiss, alpha_chain_poiss,
                                            theta_chain_poiss, M = hor, n = L)

est_ci_pred_poiss <- matrix(NA, nrow = hor, ncol = 3)
# first column = lower bound
# second columns = medians
# third columns = upper bound
for (m in 1:hor){
  est_ci_pred_poiss[m,c(1,3)] <- quantile(kmn_chain_poiss[m,], probs = c(0.025,0.975))
  est_ci_pred_poiss[m,2] <- mean(kmn_chain_poiss[m,]) 
}
est_ci_pred_poiss <- list("means" = est_ci_pred_poiss[,2],
                          "lbs" = est_ci_pred_poiss[,1],
                          "ubs" = est_ci_pred_poiss[,3])

list_kmn_pred_poiss[[labels_comb]] <- est_ci_pred_poiss

# Negative Binomial 
kmn_chain_negbin <- generate_Kmn_chain_negbin(nstar_nb, p_nb,
                                              alpha_chain_negbin, theta_chain_negbin,
                                              M = hor, n = L, Kn)

est_ci_pred_negbin <- matrix(NA, nrow = hor, ncol = 3)
# first column = lower bound
# second columns = medians
# third columns = upper bound
for (m in 1:hor){
  est_ci_pred_negbin[m,c(1,3)] <- quantile(kmn_chain_negbin[m,], probs = c(0.025,0.975))
  est_ci_pred_negbin[m,2] <- mean(kmn_chain_negbin[m,]) 
}
est_ci_pred_negbin <- list("means" = est_ci_pred_negbin[,2],
                           "lbs" = est_ci_pred_negbin[,1],
                           "ubs" = est_ci_pred_negbin[,3])

list_kmn_pred_negbin[[labels_comb]] <- est_ci_pred_negbin


# # IBP + Gamma
# kmn_chain_ibp <- generate_Kmn_chain_gamma_ibp(a_chain_ibp, b_chain_ibp, alpha_chain_ibp,
#                                               theta_chain_ibp, M = hor, n = L, Kn)
# 
# est_ci_pred_ibp <- matrix(NA, nrow = hor, ncol = 3)
# # first column = lower bound
# # second columns = medians
# # third columns = upper bound
# for (m in 1:hor){
#   est_ci_pred_ibp[m,c(1,3)] <- quantile(kmn_chain_ibp[m,], probs = c(0.025,0.975))
#   est_ci_pred_ibp[m,2] <- mean(kmn_chain_ibp[m,]) 
#   
# }
# est_ci_pred_ibp <- list("means" = est_ci_pred_ibp[,2],
#                         "lbs" = est_ci_pred_ibp[,1],
#                         "ubs" = est_ci_pred_ibp[,3])
# 
# list_kmn_pred_ibp[[labels_comb]] <- est_ci_pred_ibp

# Poisson
saveRDS(list_kmn_pred_poiss , "mazziotta2016_application/fungi/mazz_fungi_op_ci_pred_poiss.rds")
# Negative Binomial
saveRDS(list_kmn_pred_negbin, "mazziotta2016_application/fungi/mazz_fungi_op_ci_pred_negbin.rds")
# Gamma IBP
#saveRDS(list_kmn_pred_ibp, "mazziotta2016_application/fungi/mazz_fungi_op_ci_pred_ibp.rds")
