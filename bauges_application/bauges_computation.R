##############################################################
############ REAL DATA APPLICATION: Bauges ###################
#############################################################
rm(list=ls())

library(ProductFormFA)

bauges_mat <- read.csv('Bauges_results/Bauges_data.csv',header=T, row.names = 1)

freq <- colSums(bauges_mat)
# to check that all species in the dataset are actually present 
sum(freq ==0) # must be =0

# number of sites
L = nrow(bauges_mat)

bauges_mat <- bauges_mat[sample.int(L, size = L, replace = F),]

bauges_list <- create_features_list(bauges_mat)
plot_trajectory(bauges_list)

data_mat <- bauges_mat
data_list <- bauges_list

############ 1) Set parameters for the 3 models ###############

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
a_alpha_bb <- 1
b_alpha_bb <- 0.1
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


########## 1.2) Set parameters for the Gamma IBP and SB-SP

# Set prior hyperparameters for Gamma IBP
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


# Set prior hyperparameters for SB-SP
p_sp <- 0.08
print(paste0("E(c) = ", (1-p_sp)/ p_sp))
print(paste0("Var(c) = ", (1-p_sp)/ (p_sp^2)) )
r_sp <- 0.1
t_sp <- 0.01
print(paste0("E(beta) = ", r_sp/ t_sp))
print(paste0("Var(beta) = ", r_sp/ (t_sp^2)))
a_alpha_sp <- 2
b_alpha_sp <- 2
print(paste0("E(alpha) = ", a_alpha_sp/ (a_alpha_sp + b_alpha_sp)))
print(paste0("Var(alpha) = ", a_alpha_sp*b_alpha_sp/ (a_alpha_sp + b_alpha_sp)^2 /(a_alpha_sp+b_alpha_sp+1)))


# Set initial values for the parameters of Gamma IBP
a_0_ibp <- 5
b_0_ibp <- 1
alpha_0_ibp <- 0.5
s_0_ibp <- 15

# Set initial values for the parameters of SB-SP
c_0_sp <- 10
beta_0_sp <- 10
alpha_0_sp <- 0.5


########## Set MCMC parameters (common to all 3 models): 
# SB-SP has same chain settings than Gamma IBP

S_poiss <- S_negbin <- S_ibp <-  6*10^4
n_burnin_poiss <- n_burnin_negbin <- n_burnin_ibp<-  10^4
thin_poiss <- thin_negbin <- thin_ibp <- 2
seed <- 12345
number_saved_iterations_poiss <- (S_poiss - n_burnin_poiss)/thin_poiss
number_saved_iterations_negbin <- (S_negbin - n_burnin_negbin)/thin_negbin
number_saved_iterations_ibp <- (S_ibp - n_burnin_ibp)/thin_ibp

###########################################################
########## (A) Training/test approach #####################
###########################################################

Ns <- c(30, 50, 100)

###### 2) Run the algorithms ###############
gg_ntilde_poiss <- data.frame(matrix(nrow = number_saved_iterations_poiss, ncol = length(Ns)))
colnames(gg_ntilde_poiss) <- paste("N", Ns, sep = ".")
gg_ntilde_negbin <- data.frame(matrix(nrow = number_saved_iterations_negbin, ncol = length(Ns)))
colnames(gg_ntilde_negbin) <- paste("N", Ns, sep = ".")

list_kn_rarefaction_poiss <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_poiss) <- paste("N", Ns, sep = ".")
list_kn_rarefaction_negbin <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_negbin) <- paste("N", Ns, sep = ".")
list_kn_rarefaction_ibp <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_ibp) <- paste("N", Ns, sep = ".")
list_kn_rarefaction_sp <- vector(mode="list", length = length(Ns))
names(list_kn_rarefaction_sp) <- paste("N", Ns, sep = ".")

list_kmn_pred_test_poiss <- vector(mode="list", length = length(Ns))
names(list_kmn_pred_test_poiss) <- paste("N", Ns, sep = ".")
list_kmn_pred_test_negbin <- vector(mode="list", length = length(Ns))
names(list_kmn_pred_test_negbin) <- paste("N", Ns, sep = ".")
list_kmn_pred_test_ibp <- vector(mode="list", length = length(Ns))
names(list_kmn_pred_test_ibp) <- paste("N", Ns, sep = ".")
list_kmn_pred_test_sp <- vector(mode="list", length = length(Ns))
names(list_kmn_pred_test_sp) <- paste("N", Ns, sep = ".")

params_poiss <- data.frame(matrix(nrow = number_saved_iterations_poiss, ncol = 3*length(Ns)))
colnames(params_poiss) <- paste(c("lambda", "alpha", "theta"), rep(Ns, each = 3), sep = ".")
params_negbin <- data.frame(matrix(nrow = number_saved_iterations_negbin, ncol = 4*length(Ns)))
colnames(params_negbin) <- paste(c("nstar", "p", "alpha", "theta"), rep(Ns, each = 4), sep = ".")
params_ibp <- data.frame(matrix(nrow = number_saved_iterations_ibp, ncol = 4*length(Ns)))
colnames(params_ibp) <- paste(c("a", "b", "alpha", "theta"), rep(Ns, each = 4), sep = ".")
params_sp <- data.frame(matrix(nrow = number_saved_iterations_ibp, ncol = 3*length(Ns)))
colnames(params_sp) <- paste(c("c", "beta", "alpha"), rep(Ns, each = 3), sep = ".")

for (j in 1:length(Ns)){
  N <- Ns[j]
  M <- L - N
  
  train_mat <- data_mat[1:N,]
  test_mat <- data_mat[(N+1):L, ]
  # convert the binary matrix into list of features
  train_list <- create_features_list(train_mat)
  test_list <- create_features_list(test_mat)
  
  ################# 3) Run BB + Poisson ###########
  
  #Set tau for MALA
  tau_poiss <- 0.001
  
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
  tau_nb <- 0.008
  
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
  tau_ibp <- 0.005
  
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
  
  ########### 5.1) Run SB-SP #################
  
  # Set tau for MALA
  tau_sp <- 0.005
  
  output_sp <- gibbs_sampler_sb_sp(Z = train_mat,
                                   c_0_sp, beta_0_sp, alpha_0_sp,
                                   p_sp, r_sp, t_sp, a_alpha_sp, b_alpha_sp,
                                   tau_sp, fixed = c(F,F,F),
                                   S_ibp, n_burnin_ibp, thin_ibp, seed)
  
  n_saved_iter_sp <- length(output_sp$c_vec)
  c_chain_sp <- output_sp$c_vec
  beta_chain_sp <- output_sp$beta_vec
  alpha_chain_sp <- output_sp$alpha_vec
  
  params_sp[[paste0("c.",N)]] <- c_chain_sp
  params_sp[[paste0("beta.",N)]] <- beta_chain_sp
  params_sp[[paste0("alpha.",N)]] <- alpha_chain_sp
  
  ####### 6) Estimate limit distributions (Poiss/NB) ################
  Kn = ncol(train_mat[,colSums(train_mat) > 0])
  
  ntilde_chain_poiss <- generate_Ntilde_chain_poiss(lambda_chain_poiss, alpha_chain_poiss,
                                                    theta_chain_poiss, n = N, Kn)
  
  ntilde_chain_negbin <- generate_Ntilde_chain_negbin(nstar_chain_negbin, p_chain_negbin,
                                                      alpha_chain_negbin,
                                                      theta_chain_negbin, n = N, Kn)
  
  gg_ntilde_poiss[[paste0("N.",N)]] <- ntilde_chain_poiss
  gg_ntilde_negbin[[paste0("N.",N)]] <- ntilde_chain_negbin
  
  ######## 7) Extrapolation in the test set (Poiss/NB/Gamma/SP) ##############
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
  
  # SB-SP
  kmn_chain_sp <- generate_Kmn_chain_gamma_ibp(a_chain = c_chain_sp + 1,
                                               b_chain = beta_chain_sp*(1-alpha_chain_sp)/alpha_chain_sp, 
                                               alpha_chain = alpha_chain_sp,
                                               theta_chain = 1 - alpha_chain_sp,
                                               M = M, n = N, Kn)
  
  est_ci_pred_sp <- matrix(NA, nrow = M, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:M){
    est_ci_pred_sp[m,] <- quantile(kmn_chain_sp[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_pred_sp <- list("medians" = est_ci_pred_sp[,2],
                         "lbs" = est_ci_pred_sp[,1],
                         "ubs" = est_ci_pred_sp[,3])
  
  list_kmn_pred_test_sp[[paste0("N.",N)]] <- est_ci_pred_sp
  
  ############ 7.b) Rarefaction curve ################
  # Poisson
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
  
  # Neg-Bin
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
  
  # SB-SP
  kn_chain_sp <- generate_Kmn_chain_gamma_ibp(a_chain = c_chain_sp + 1,
                                              b_chain = beta_chain_sp*(1-alpha_chain_sp)/alpha_chain_sp, 
                                              alpha_chain = alpha_chain_sp,
                                              theta_chain = 1 - alpha_chain_sp,
                                              M = N, n = 0)
  
  est_ci_sp <- matrix(NA, nrow = N, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:N){
    est_ci_sp[m,] <- quantile(kn_chain_sp[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_sp <- list("medians" = est_ci_sp[,2],
                    "lbs" = est_ci_sp[,1],
                    "ubs" = est_ci_sp[,3])
  
  list_kn_rarefaction_sp[[paste0("N.",N)]] <- est_ci_sp
  
}


############# 8) Save results:  MCMC convergence ####################

save(params_poiss, file = "bauges_results/bauges_params_poiss.Rda")
save(params_negbin, file = "bauges_results/bauges_params_negbin.Rda")
save(params_ibp, file = "bauges_results/bauges_params_ibp.Rda")
save(params_sp, file = "bauges_results/bauges_params_sp.Rda")

############ 9) Save results: samples from limiting distributions (Poiss/NB) ##############
# Poisson
save(gg_ntilde_poiss, file = "bauges_results/bauges_ntilde_poiss.Rda")
# Negative Binomial
save(gg_ntilde_negbin, file = "bauges_results/bauges_ntilde_negbin.Rda")

######### 10) Save results: CI for extrapolation (Poiss/NB/Gamma/SP) ################
# Poisson
saveRDS(list_kmn_pred_test_poiss, "bauges_results/bauges_ci_poiss.rds")
# Negative Binomial
saveRDS(list_kmn_pred_test_negbin, "bauges_results/bauges_ci_negbin.rds")
# Gamma IBP
saveRDS(list_kmn_pred_test_ibp, "bauges_results/bauges_ci_ibp.rds")
# SB-SP
saveRDS(list_kmn_pred_test_sp, "bauges_results/bauges_ci_sp.rds")

######### 10.b) Save results: CI for rarefaction (Poiss/NB/Gamma/SP) ################
# Poisson
saveRDS(list_kn_rarefaction_poiss, "bauges_results/bauges_ci_insample_poiss.rds")
# Negative Binomial
saveRDS(list_kn_rarefaction_negbin, "bauges_results/bauges_ci_insample_negbin.rds")
# Gamma IBP
saveRDS(list_kn_rarefaction_ibp, "bauges_results/bauges_ci_insample_ibp.rds")
# SB-SP
saveRDS(list_kn_rarefaction_sp, "bauges_results/bauges_ci_insample_sp.rds")

######## 11) Save the data ##############################
saveRDS(data_mat, "bauges_results/bauges_data_mat.rds")




###########################################################
########## (B) Only prediction - NOT SO INTERESTING HERE SINCE ######
######### THE TEST SET IS SATURATING ALREADY #####################
###### LOOK JUST AT TRAINING/TEST APPROACH ######################
###########################################################

###### 2) Run the algorithms ###############
gg_ntilde_poiss <- data.frame(matrix(nrow = number_saved_iterations_poiss, ncol = 1))
colnames(gg_ntilde_poiss) <- paste("N", L, sep = ".")
gg_ntilde_negbin <- data.frame(matrix(nrow = number_saved_iterations_negbin, ncol = 1))
colnames(gg_ntilde_negbin) <- paste("N", L, sep = ".")

params_poiss <- data.frame(matrix(nrow = number_saved_iterations_poiss, ncol = 3))
colnames(params_poiss) <- c("lambda", "alpha", "theta")
params_negbin <- data.frame(matrix(nrow = number_saved_iterations_negbin, ncol = 4))
colnames(params_negbin) <- c("nstar", "p", "alpha", "theta")
params_ibp <- data.frame(matrix(nrow = number_saved_iterations_ibp, ncol = 4))
colnames(params_ibp) <- c("a", "b", "alpha", "theta")
params_sp <- data.frame(matrix(nrow = number_saved_iterations_ibp, ncol = 3))
colnames(params_sp) <- c("c", "beta", "alpha")


################# 3) Run BB + Poisson ###########

#Set tau for MALA
tau_poiss <- 0.007

output_poiss <- gibbs_sampler_poiss(Z = data_mat,
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

params_poiss[["lambda"]] <- lambda_chain_poiss
params_poiss[["alpha"]] <- alpha_chain_poiss
params_poiss[["theta"]] <- theta_chain_poiss


####### 4) Run BB + Negative-Binomial ################

# Set tau for MALA
tau_nb <- 0.01

output_negbin <- gibbs_sampler_negbin_geometric(Z = data_mat,
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

params_negbin[["nstar"]] <- nstar_chain_negbin
params_negbin[["p"]] <- p_chain_negbin
params_negbin[["alpha"]] <- alpha_chain_negbin
params_negbin[["theta"]] <- theta_chain_negbin

############ 5) Run IBP + Gamma #################

# Set tau for MALA
tau_ibp <- 0.008

output_ibp <- gibbs_sampler_gamma_ibp(Z = data_mat,
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

params_ibp[["a"]] <- a_chain_ibp
params_ibp[["b"]] <- b_chain_ibp
params_ibp[["alpha"]] <- alpha_chain_ibp
params_ibp[["theta"]] <- theta_chain_ibp

########### 5.1) Run SB-SP #################

# Set tau for MALA
tau_sp <- 0.02

output_sp <- gibbs_sampler_sb_sp(Z = data_mat,
                                 c_0_sp, beta_0_sp, alpha_0_sp,
                                 p_sp, r_sp, t_sp, a_alpha_sp, b_alpha_sp,
                                 tau_sp, fixed = c(F,F,F),
                                 S_ibp, n_burnin_ibp, thin_ibp, seed)

n_saved_iter_sp <- length(output_sp$c_vec)
c_chain_sp <- output_sp$c_vec
beta_chain_sp <- output_sp$beta_vec
alpha_chain_sp <- output_sp$alpha_vec

params_sp[["c"]] <- c_chain_sp
params_sp[["beta"]] <- beta_chain_sp
params_sp[["alpha"]] <- alpha_chain_sp


####### 6) Estimate limit distributions (Poiss/NB) ################
Kn = ncol(data_mat[,colSums(data_mat) > 0])

ntilde_chain_poiss <- generate_Ntilde_chain_poiss(lambda_chain_poiss, alpha_chain_poiss,
                                                  theta_chain_poiss, n = L, Kn)

ntilde_chain_negbin <- generate_Ntilde_chain_negbin(nstar_chain_negbin, p_chain_negbin,
                                                    alpha_chain_negbin,
                                                    theta_chain_negbin, n = L, Kn)

gg_ntilde_poiss[[paste0("N.",L)]] <- ntilde_chain_poiss
gg_ntilde_negbin[[paste0("N.",L)]] <- ntilde_chain_negbin


############ 7) Rarefaction curve ################

# Poisson
kn_chain_poiss <- generate_Kmn_chain_poiss(lambda_chain_poiss, alpha_chain_poiss,
                                           theta_chain_poiss, M = L, n=0)

est_ci_poiss <- matrix(NA, nrow = L, ncol = 3)
# first column = lower bound
# second columns = medians
# third columns = upper bound
for (m in 1:L){
  est_ci_poiss[m,] <- quantile(kn_chain_poiss[m,], probs = c(0.025,0.5,0.975))
}
est_ci_poiss <- list("medians" = est_ci_poiss[,2],
                     "lbs" = est_ci_poiss[,1],
                     "ubs" = est_ci_poiss[,3])

kn_rarefaction_poiss <- est_ci_poiss

# Neg-Bin
kn_chain_negbin <- generate_Kmn_chain_negbin(nstar_chain_negbin, p_chain_negbin,
                                             alpha_chain_negbin, theta_chain_negbin,
                                             M = L, n=0)

est_ci_negbin <- matrix(NA, nrow = L, ncol = 3)
# first column = lower bound
# second columns = medians
# third columns = upper bound
for (m in 1:L){
  est_ci_negbin[m,] <- quantile(kn_chain_negbin[m,], probs = c(0.025,0.5,0.975))
}
est_ci_negbin <- list("medians" = est_ci_negbin[,2],
                      "lbs" = est_ci_negbin[,1],
                      "ubs" = est_ci_negbin[,3])

kn_rarefaction_negbin <- est_ci_negbin

# IBP
kn_chain_ibp <- generate_Kmn_chain_gamma_ibp(a_chain_ibp, b_chain_ibp, alpha_chain_ibp, 
                                             theta_chain_ibp, M = L, n=0)

est_ci_ibp <- matrix(NA, nrow = L, ncol = 3)
# first column = lower bound
# second columns = medians
# third columns = upper bound
for (m in 1:L){
  est_ci_ibp[m,] <- quantile(kn_chain_ibp[m,], probs = c(0.025,0.5,0.975))
}
est_ci_ibp <- list("medians" = est_ci_ibp[,2],
                   "lbs" = est_ci_ibp[,1],
                   "ubs" = est_ci_ibp[,3])

kn_rarefaction_ibp <- est_ci_ibp

# SB-SP
kn_chain_sp <- generate_Kmn_chain_gamma_ibp(a_chain = c_chain_sp + 1,
                                            b_chain = beta_chain_sp*(1-alpha_chain_sp)/alpha_chain_sp, 
                                            alpha_chain = alpha_chain_sp,
                                            theta_chain = 1 - alpha_chain_sp,
                                            M = L, n = 0)

est_ci_sp <- matrix(NA, nrow = L, ncol = 3)
# first column = lower bound
# second columns = medians
# third columns = upper bound
for (m in 1:L){
  est_ci_sp[m,] <- quantile(kn_chain_sp[m,], probs = c(0.025,0.5,0.975))
}
est_ci_sp <- list("medians" = est_ci_sp[,2],
                  "lbs" = est_ci_sp[,1],
                  "ubs" = est_ci_sp[,3])

kn_rarefaction_sp <- est_ci_sp





############# 8) Save results:  MCMC convergence ####################

save(params_poiss, file = "bauges_results/bauges_op_params_poiss.Rda")
save(params_negbin, file = "bauges_results/bauges_op_params_negbin.Rda")
save(params_ibp, file = "bauges_results/bauges_op_params_ibp.Rda")
save(params_sp, file = "bauges_results/bauges_op_params_sp.Rda")

############ 9) Save results: samples from limiting distributions (Poiss/NB) ##############
# Poisson
save(gg_ntilde_poiss, file = "bauges_results/bauges_op_ntilde_poiss.Rda")
# Negative Binomial
save(gg_ntilde_negbin, file = "bauges_results/bauges_op_ntilde_negbin.Rda")

######### 10.b) Save results: CI for rarefaction (Poiss/NB/Gamma/SP) ################
# Poisson
saveRDS(kn_rarefaction_poiss, "bauges_results/bauges_op_ci_insample_poiss.rds")
# Negative Binomial
saveRDS(kn_rarefaction_negbin, "bauges_results/bauges_op_ci_insample_negbin.rds")
# Gamma IBP
saveRDS(kn_rarefaction_ibp, "bauges_results/bauges_op_ci_insample_ibp.rds")
# SB-SP
saveRDS(kn_rarefaction_sp, "bauges_results/bauges_op_ci_insample_sp.rds")

######## 11) Save the data ##############################
saveRDS(data_mat, "bauges_results/bauges_data_mat.rds")
