####
##### MODEL 5: the Zipf-Mandelbrot model ############################
####

rm(list=ls())
library(tidyverse)

# set number of individuals
Ns <- c(20, 40, 60, 80, 150) 
L <- 500

# maximum number of features
H <- 500
# choose c such that pi's <= 0.5, where pi_i = c/(i+5), i=1,...,H
c <- 3
# define pi's
pis <- c/((1:H) + 5)

###### 1) Set parameters for the 3 models ###############

########### 1.1) Set parameters for BB (Poisson and Neg-Bin)

# Set prior hyperparameters specific for poisson
a_l_poiss <- 50
b_l_poiss <- 0.1
print(paste0("E(lambda) = ", a_l_poiss/ b_l_poiss))
print(paste0("Var(lambda) = ", a_l_poiss/ (b_l_poiss^2)))

# Set prior hyperparameters specific for NB
q_star_nb <- 0.01
print(paste0("E(nstar) = ", 1/ q_star_nb))
print(paste0("Var(nstar) = ", (1-q_star_nb)/ (q_star_nb^2)) )
alpha_p_nb <- 0.1
beta_p_nb <- 0.1
print(paste0("E(p) = ", alpha_p_nb/(alpha_p_nb + beta_p_nb)))
print(paste0("Var(p) = ",(alpha_p_nb*beta_p_nb)/((alpha_p_nb + beta_p_nb)^2 *(alpha_p_nb+beta_p_nb+1)) ))

# Set other hyperparameters
a_alpha_bb <- 5
b_alpha_bb <- 0.5
print(paste0("E(alpha_bar) = ", a_alpha_bb/ b_alpha_bb))
print(paste0("Var(alpha_bar) = ", a_alpha_bb/ (b_alpha_bb^2)))
a_s_bb <- 2
b_s_bb <- 0.2
print(paste0("E(s) = ", a_s_bb/ b_s_bb))
print(paste0("Var(s) = ", a_s_bb/ (b_s_bb^2)))

# Set initial values for poisson 
lambda_0_poiss <- 100
# Set initial values for NB
nstar_0_nb <- 100
p_0_nb <- 0.2
# Set initial values for other parameters
alpha_bar_0_bb <- 1
s_0_bb <- 1


########## 1.2) Set parameters for the Gamma IBP

# Set prior hyperparameters
p_ibp <- 0.05
print(paste0("E(a) = ", 1/ p_ibp))
print(paste0("Var(a) = ", (1-p_ibp)/ (p_ibp^2)) )
r_ibp <- 1
t_ibp <- 0.1
print(paste0("E(b) = ", r_ibp/ t_ibp))
print(paste0("Var(b) = ", r_ibp/ (t_ibp^2)))
a_alpha_ibp <- 2
b_alpha_ibp <- 2
print(paste0("E(alpha) = ", a_alpha_ibp/ (a_alpha_ibp + b_alpha_ibp)))
print(paste0("Var(alpha) = ", a_alpha_ibp*b_alpha_ibp/ (a_alpha_ibp + b_alpha_ibp)^2 /(a_alpha_ibp+b_alpha_ibp+1)))
a_s_ibp <- 2
b_s_ibp <- 0.2
print(paste0("E(s) = ", a_s_ibp/ b_s_ibp))
print(paste0("Var(s) = ", a_s_ibp/ (b_s_ibp^2)))

print(paste0("E(theta) = ", a_s_ibp/ b_s_ibp -  a_alpha_ibp/ (a_alpha_ibp + b_alpha_ibp)))
print(paste0("Var(theta) = ", a_s_ibp/ (b_s_ibp^2) + 
               a_alpha_ibp*b_alpha_ibp/ (a_alpha_ibp + b_alpha_ibp)^2 /(a_alpha_ibp+b_alpha_ibp+1)))


# Set initial values for the parameters of Gamma IBP
a_0_ibp <- 5
b_0_ibp <- 1
alpha_0_ibp <- 0.5
s_0_ibp <- 15


##### Single dataset -> Ntilde (Poiss/NB) and extrapolation (Poiss/NB/Gamma) #####
seed = 123456
set.seed(seed)
data_mat <- matrix(rbinom(L*H, size = 1, prob = rep(pis, L)), 
                   nrow = L, ncol = H, byrow = T )
data_list <- create_features_list(data_mat)
#plot_trajectory(data_list)

num_feat <- vector(length = length(Ns))
for (j in 1:length(Ns)){
  N <- Ns[j]
  num_feat[j] <- sum(colSums(data_mat[1:N, ]) > 0)
}
#print("N. of observed features in the sample: ")
#print(num_feat)

########## Set MCMC parameters (common to all 3 models)

S_poiss <- S_negbin <- S_ibp <- 6*10^4
n_burnin_poiss <- n_burnin_negbin <- n_burnin_ibp<- 10^4
thin_poiss <- thin_negbin <- thin_ibp <- 2
seed <- 1234
number_saved_iterations_poiss <- (S_poiss - n_burnin_poiss)/thin_poiss 
number_saved_iterations_negbin <- (S_negbin - n_burnin_negbin)/thin_negbin 
number_saved_iterations_ibp <- (S_ibp - n_burnin_ibp)/thin_ibp 

###### 2) Run the algorithms ###############
gg_ntilde_poiss <- data.frame(matrix(nrow = number_saved_iterations_poiss, ncol = length(Ns)))
colnames(gg_ntilde_poiss) <- paste("N", Ns, sep = ".")
gg_ntilde_negbin <- data.frame(matrix(nrow = number_saved_iterations_negbin, ncol = length(Ns)))
colnames(gg_ntilde_negbin) <- paste("N", Ns, sep = ".")

list_kmn_pred_test_poiss <- vector(mode="list", length = length(Ns))
names(list_kmn_pred_test_poiss) <- paste("N", Ns, sep = ".")
list_kmn_pred_test_negbin <- vector(mode="list", length = length(Ns))
names(list_kmn_pred_test_negbin) <- paste("N", Ns, sep = ".")
list_kmn_pred_test_ibp <- vector(mode="list", length = length(Ns))
names(list_kmn_pred_test_ibp) <- paste("N", Ns, sep = ".")

params_poiss <- data.frame(matrix(nrow = number_saved_iterations_poiss, ncol = 3*length(Ns)))
colnames(params_poiss) <- paste(c("lambda", "alpha", "theta"), rep(Ns, each = 3), sep = ".")
params_negbin <- data.frame(matrix(nrow = number_saved_iterations_negbin, ncol = 4*length(Ns)))
colnames(params_negbin) <- paste(c("nstar", "p", "alpha", "theta"), rep(Ns, each = 4), sep = ".")
params_ibp <- data.frame(matrix(nrow = number_saved_iterations_ibp, ncol = 4*length(Ns)))
colnames(params_ibp) <- paste(c("a", "b", "alpha", "theta"), rep(Ns, each = 4), sep = ".")

for (j in 1:length(Ns)){
  N <- Ns[j]
  M <- L - N
  
  train_mat <- data_mat[1:N,]
  test_mat <- data_mat[(N+1):L, ]
  # convert the binary matrix into list of features
  train_list <- create_features_list(train_mat)
  test_list <- create_features_list(test_mat)
  
  ################# 3) Run BB + Poisson ###########
  
  # Set tau for MALA
  tau_poiss <- 0.0005
  
  output_poiss <- gibbs_sampler_poiss(Z = train_mat,
                                      lambda_0_poiss, alpha_bar_0_bb, s_0_bb,
                                      a_l_poiss, b_l_poiss, a_alpha_bb, b_alpha_bb, a_s_bb, b_s_bb,
                                      tau_poiss,
                                      S_poiss, n_burnin_poiss, thin_poiss, seed)
  
  n_saved_iter_poiss <- length(output_poiss$lambda_vec)
  lambda_chain_poiss <- output_poiss$lambda_vec
  s_chain_poiss <- output_poiss$s_vec
  alpha_bar_chain_poiss <- output_poiss$alpha_bar_vec
  alpha_chain_poiss <- - alpha_bar_chain_poiss
  theta_chain_poiss <- s_chain_poiss + alpha_bar_chain_poiss
  
  params_poiss[[paste0("lambda.",N)]] <- lambda_chain_poiss
  params_poiss[[paste0("alpha.",N)]] <- alpha_chain_poiss
  params_poiss[[paste0("theta.",N)]] <- theta_chain_poiss
  
  
  ####### 4) Run BB + Negative-Binomial ################
  
  # Set tau for MALA
  tau_nb <- 0.002
  
  output_negbin <- gibbs_sampler_negbin_geometric(Z = train_mat,
                                                  nstar_0_nb, p_0_nb,  s_0_bb, alpha_bar_0_bb,
                                                  q_star_nb, alpha_p_nb, beta_p_nb, 
                                                  a_alpha_bb, b_alpha_bb, a_s_bb, b_s_bb,
                                                  tau_nb, fixed = c(F,F,F,F),
                                                  S_negbin, n_burnin_negbin, thin_negbin, seed)
  
  n_saved_iter_negbin <- length(output_negbin$nstar_vec)
  nstar_chain_negbin <- output_negbin$nstar_vec
  p_chain_negbin <- output_negbin$p_vec
  s_chain_negbin <- output_negbin$s_vec
  alpha_bar_chain_negbin <- output_negbin$alpha_bar_vec
  alpha_chain_negbin <- - alpha_bar_chain_negbin
  theta_chain_negbin <- s_chain_negbin + alpha_bar_chain_negbin
  
  params_negbin[[paste0("nstar.",N)]] <- nstar_chain_negbin
  params_negbin[[paste0("p.",N)]] <- p_chain_negbin
  params_negbin[[paste0("alpha.",N)]] <- alpha_chain_negbin
  params_negbin[[paste0("theta.",N)]] <- theta_chain_negbin
  
  ############ 5) Run IBP + Gamma #################
  
  # Set tau for MALA
  tau_ibp <- 0.002
  
  output_ibp <- gibbs_sampler_gamma_ibp(Z = train_mat,
                                        a_0_ibp, b_0_ibp, s_0_ibp, alpha_0_ibp,
                                        p_ibp, r_ibp, t_ibp, a_alpha_ibp, b_alpha_ibp, a_s_ibp, b_s_ibp,
                                        tau_ibp, fixed = c(F,F,F,F),
                                        S_ibp, n_burnin_ibp, thin_ibp, seed)
  
  n_saved_iter_ibp <- length(output_ibp$a_vec)
  a_chain_ibp <- output_ibp$a_vec
  b_chain_ibp <- output_ibp$b_vec
  s_chain_ibp <- output_ibp$s_vec
  alpha_chain_ibp <- output_ibp$alpha_vec
  theta_chain_ibp <- s_chain_ibp - alpha_chain_ibp
  
  params_ibp[[paste0("a.",N)]] <- a_chain_ibp
  params_ibp[[paste0("b.",N)]] <- b_chain_ibp
  params_ibp[[paste0("alpha.",N)]] <- alpha_chain_ibp
  params_ibp[[paste0("theta.",N)]] <- theta_chain_ibp
  
  ####### 6) Estimate limit distributions (Poiss/NB) ################
  Kn = ncol(train_mat[,colSums(train_mat) > 0])
  
  ntilde_chain_poiss <- generate_Ntilde_chain_poiss(lambda_chain_poiss, alpha_chain_poiss, 
                                                    theta_chain_poiss, n = N, Kn)
  #ntilde_chain_poiss <- as.data.frame(ntilde_chain_poiss)
  
  ntilde_chain_negbin <- generate_Ntilde_chain_negbin(nstar_chain_negbin, p_chain_negbin,
                                                      alpha_chain_negbin,
                                                      theta_chain_negbin, n = N, Kn)
  #ntilde_chain_negbin <- as.data.frame(ntilde_chain_negbin)
  
  # gg_ntilde_poiss[[j]] <- ggplot(ntilde_chain_poiss, aes(x=ntilde_chain_poiss)) + 
  #   geom_density() + ggtitle(paste0("Poisson, N = ", N)) 
  # gg_ntilde_negbin[[j]] <- ggplot(ntilde_chain_negbin, aes(x=ntilde_chain_negbin)) + 
  #   geom_density() + ggtitle(paste0("Negbin, N = ", N))
  
  gg_ntilde_poiss[[paste0("N.",N)]] <- ntilde_chain_poiss
  gg_ntilde_negbin[[paste0("N.",N)]] <- ntilde_chain_negbin
  
  ######## 7) Extrapolation in the test set (Poiss/NB/Gamma) ##############
  # Poisson
  kmn_chain_poiss <- generate_Kmn_chain_poiss(lambda_chain_poiss, alpha_chain_poiss,
                                              theta_chain_poiss, M = M, n = N)
  
  est_ci_pred_poiss <- matrix(NA, nrow = M, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:M){
    est_ci_pred_poiss[m,] <- quantile(kmn_chain_poiss[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_pred_poiss <- list("medians" = est_ci_pred_poiss[,2],
                            "lbs" = est_ci_pred_poiss[,1],
                            "ubs" = est_ci_pred_poiss[,3])
  
  list_kmn_pred_test_poiss[[paste0("N.",N)]] <- est_ci_pred_poiss
  
  # Negative Binomial
  kmn_chain_negbin <- generate_Kmn_chain_negbin(nstar_chain_negbin, p_chain_negbin,
                                                alpha_chain_negbin, theta_chain_negbin,
                                                M = M, n = N, Kn)
  
  est_ci_pred_negbin <- matrix(NA, nrow = M, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:M){
    est_ci_pred_negbin[m,] <- quantile(kmn_chain_negbin[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_pred_negbin <- list("medians" = est_ci_pred_negbin[,2],
                             "lbs" = est_ci_pred_negbin[,1],
                             "ubs" = est_ci_pred_negbin[,3])
  
  list_kmn_pred_test_negbin[[paste0("N.",N)]] <- est_ci_pred_negbin
  
  
  # IBP + Gamma
  kmn_chain_ibp <- generate_Kmn_chain_gamma_ibp(a_chain_ibp, b_chain_ibp, alpha_chain_ibp, 
                                                theta_chain_ibp, M = M, n = N, Kn)
  
  est_ci_pred_ibp <- matrix(NA, nrow = M, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:M){
    est_ci_pred_ibp[m,] <- quantile(kmn_chain_ibp[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_pred_ibp <- list("medians" = est_ci_pred_ibp[,2],
                          "lbs" = est_ci_pred_ibp[,1],
                          "ubs" = est_ci_pred_ibp[,3])
  
  list_kmn_pred_test_ibp[[paste0("N.",N)]] <- est_ci_pred_ibp
  
}


############# 8) Save results:  MCMC convergence ####################

save(params_poiss, file = "chao_model_simulation/m5/m5_params_poiss.Rda")
save(params_negbin, file = "chao_model_simulation/m5/m5_params_negbin.Rda")
save(params_ibp, file = "chao_model_simulation/m5/m5_params_ibp.Rda")

############ 9) Save results: samples from limiting distributions (Poiss/NB) ##############
# Poisson
save(gg_ntilde_poiss, file = "chao_model_simulation/m5/m5_ntilde_poiss.Rda")
# Negative Binomial
save(gg_ntilde_negbin, file = "chao_model_simulation/m5/m5_ntilde_negbin.Rda")

######### 10) Save results: CI for extrapolation (Poiss/NB/Gamma) ################
# Poisson
saveRDS(list_kmn_pred_test_poiss, "chao_model_simulation/m5/m5_ci_poiss.rds")
# Negative Binomial
saveRDS(list_kmn_pred_test_negbin, "chao_model_simulation/m5/m5_ci_negbin.rds")
# Gamma IBP
saveRDS(list_kmn_pred_test_ibp, "chao_model_simulation/m5/m5_ci_ibp.rds")


######## 11) Save the data ##############################
saveRDS(data_mat, "chao_model_simulation/m5/m5_data_mat.rds")




##### Accuracy on multiple datasets #####

# number of datasets to average over
D <- 50

# mcmc parameters
S_poiss <- S_negbin <- S_ibp <- 5*10^4
n_burnin_poiss <- n_burnin_negbin <- n_burnin_ibp<- 10^4
thin_poiss <- thin_negbin <- thin_ibp <- 2
number_saved_iterations_poiss <- (S_poiss - n_burnin_poiss)/thin_poiss 
number_saved_iterations_negbin <- (S_negbin - n_burnin_negbin)/thin_negbin 
number_saved_iterations_ibp <- (S_ibp - n_burnin_ibp)/thin_ibp 

# store objects
avg_ntilde_poiss <- data.frame(matrix(nrow = D, ncol = length(Ns)))
colnames(avg_ntilde_poiss) <- paste("N", Ns, sep = ".")
avg_ntilde_negbin <- data.frame(matrix(nrow = D, ncol = length(Ns)))
colnames(avg_ntilde_negbin) <- paste("N", Ns, sep = ".")

# number of new features observed in the test
obs_new <- data.frame(matrix(nrow = D, ncol = length(Ns)))
colnames(obs_new) <- paste("N", Ns, sep = ".")

# number of old features observed in the training
obs_train <- data.frame(matrix(nrow = D, ncol = length(Ns)))
colnames(obs_train) <- paste("N", Ns, sep = ".")

est_new_poiss <- data.frame(matrix(nrow = D, ncol = length(Ns)))
colnames(est_new_poiss) <- paste("N", Ns, sep = ".")
est_new_negbin <- data.frame(matrix(nrow = D, ncol = length(Ns)))
colnames(est_new_negbin) <- paste("N", Ns, sep = ".")
est_new_ibp <- data.frame(matrix(nrow = D, ncol = length(Ns)))
colnames(est_new_ibp) <- paste("N", Ns, sep = ".")

seed = 123456
set.seed(seed)

for (d in 1:D){
  data_mat <- matrix(rbinom(L*H, size = 1, prob = rep(pis, L)), 
                     nrow = L, ncol = H, byrow = T )
  data_list <- create_features_list(data_mat)
  
  for (j in 1:length(Ns)){
    N <- Ns[j]
    M <- L - N
    
    train_mat <- data_mat[1:N,]
    test_mat <- data_mat[(N+1):L, ]
    # convert the binary matrix into list of features
    train_list <- create_features_list(train_mat)
    test_list <- create_features_list(test_mat)
    
    ###### 3) Run BB + Poisson ###########
    
    # Set tau for MALA
    tau_poiss <- 0.0005
    
    output_poiss <- gibbs_sampler_poiss(Z = train_mat,
                                        lambda_0_poiss, alpha_bar_0_bb, s_0_bb,
                                        a_l_poiss, b_l_poiss, a_alpha_bb, b_alpha_bb, a_s_bb, b_s_bb,
                                        tau_poiss,
                                        S_poiss, n_burnin_poiss, thin_poiss, seed)
    
    n_saved_iter_poiss <- length(output_poiss$lambda_vec)
    lambda_chain_poiss <- output_poiss$lambda_vec
    s_chain_poiss <- output_poiss$s_vec
    alpha_bar_chain_poiss <- output_poiss$alpha_bar_vec
    alpha_chain_poiss <- - alpha_bar_chain_poiss
    theta_chain_poiss <- s_chain_poiss + alpha_bar_chain_poiss
    
    ##### 4) Run BB + Negative-Binomial ################
    
    # Set tau for MALA
    tau_nb <- 0.001
    
    output_negbin <- gibbs_sampler_negbin_geometric(Z = train_mat,
                                                    nstar_0_nb, p_0_nb,  s_0_bb, alpha_bar_0_bb,
                                                    q_star_nb, alpha_p_nb, beta_p_nb, 
                                                    a_alpha_bb, b_alpha_bb, a_s_bb, b_s_bb,
                                                    tau_nb, fixed = c(F,F,F,F),
                                                    S_negbin, n_burnin_negbin, thin_negbin, seed)
    
    n_saved_iter_negbin <- length(output_negbin$nstar_vec)
    nstar_chain_negbin <- output_negbin$nstar_vec
    p_chain_negbin <- output_negbin$p_vec
    s_chain_negbin <- output_negbin$s_vec
    alpha_bar_chain_negbin <- output_negbin$alpha_bar_vec
    alpha_chain_negbin <- - alpha_bar_chain_negbin
    theta_chain_negbin <- s_chain_negbin + alpha_bar_chain_negbin
    
    
    ##### 5) Run IBP + Gamma #################
    
    # Set tau for MALA
    tau_ibp <- 0.002
    
    output_ibp <- gibbs_sampler_gamma_ibp(Z = train_mat,
                                          a_0_ibp, b_0_ibp, s_0_ibp, alpha_0_ibp,
                                          p_ibp, r_ibp, t_ibp, a_alpha_ibp, b_alpha_ibp, a_s_ibp, b_s_ibp,
                                          tau_ibp, fixed = c(F,F,F,F),
                                          S_ibp, n_burnin_ibp, thin_ibp, seed)
    
    n_saved_iter_ibp <- length(output_ibp$a_vec)
    a_chain_ibp <- output_ibp$a_vec
    b_chain_ibp <- output_ibp$b_vec
    s_chain_ibp <- output_ibp$s_vec
    alpha_chain_ibp <- output_ibp$alpha_vec
    theta_chain_ibp <- s_chain_ibp - alpha_chain_ibp
    
    
    ##### 6) Estimate mean of the limit distributions (Poiss/NB) ################
    Kn = ncol(train_mat[,colSums(train_mat) > 0])
    
    ntilde_chain_poiss <- generate_Ntilde_chain_poiss(lambda_chain_poiss, alpha_chain_poiss, 
                                                      theta_chain_poiss, n = N, Kn)
    
    ntilde_chain_negbin <- generate_Ntilde_chain_negbin(nstar_chain_negbin, p_chain_negbin,
                                                        alpha_chain_negbin,
                                                        theta_chain_negbin, n = N, Kn)
    
    avg_ntilde_poiss[d, j] <- mean(ntilde_chain_poiss)
    avg_ntilde_negbin[d, j] <- mean(ntilde_chain_negbin)
    
    ##### 7) Collect quantities for accuracy (Poiss/NB/Gamma) #####
    feat_train <- unique(unlist(train_list))
    feat_test <- unique(unlist(test_list))
    # number of features observed in training
    obs_train[d,j] <- length(feat_train)
    # number of new features observed in test
    obs_new_features <- setdiff(feat_test, feat_train)
    obs_new[d,j] <- length(obs_new_features)
    
    # Poisson
    kmn_chain_poiss <- generate_Kmn_chain_poiss(lambda_chain_poiss, alpha_chain_poiss,
                                                theta_chain_poiss, M = M, n = N)
    
    est_new_poiss[d , j] <- mean(kmn_chain_poiss)
    
    # Negative Binomial
    kmn_chain_negbin <- generate_Kmn_chain_negbin(nstar_chain_negbin, p_chain_negbin,
                                                  alpha_chain_negbin, theta_chain_negbin,
                                                  M = M, n = N, Kn)
    
    est_new_negbin[d , j] <- mean(kmn_chain_negbin)
    
    # IBP + Gamma
    kmn_chain_ibp <- generate_Kmn_chain_gamma_ibp(a_chain_ibp, b_chain_ibp, alpha_chain_ibp, 
                                                  theta_chain_ibp, M = M, n = N, Kn)
    
    est_new_ibp[d , j] <- mean(kmn_chain_ibp)
    
  }
  
}


###### 8) Save results: limit distribution estimates #####
save(avg_ntilde_poiss, file = "chao_model_simulation/m5/m5_avg_ntilde_poiss.Rda")
save(avg_ntilde_negbin, file = "chao_model_simulation/m5/m5_avg_ntilde_negbin.Rda")


###### 9) Save results: quantities for accuracy #####
save(obs_train, file = "chao_model_simulation/m5/m5_obs_train.Rda")
save(obs_new, file = "chao_model_simulation/m5/m5_obs_new.Rda")
save(est_new_poiss, file = "chao_model_simulation/m5/m5_est_new_poiss.Rda")
save(est_new_negbin, file = "chao_model_simulation/m5/m5_est_new_negbin.Rda")
save(est_new_ibp, file = "chao_model_simulation/m5/m5_est_new_ibp.Rda")

