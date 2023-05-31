##############################################################
###### Simulation in FINITE-DIMENSIONAL CASE #################
##############################################################
rm(list=ls())
library(ggpubr)

#############################################################
################# SETTING 2: ZIPF (linear growth ) ########################
#############################################################

############ 1) Generate the data #####################

# set number of individuals
L <- 1000
N <- 300 # dimension of the training (used as whole sample when we don't consider test/train)
M <- L - N
# set maximum number of features in the population
H <- 10^5

# value of xi
xi = 1.2

seed = 1234
# generate data from zipf, in form of binary matrix
data_mat <- rzipf(n = L, K = H, xi = xi, seed = seed)

# convert the binary matrix into list of features
data_list <- create_features_list(data_mat)
plot_trajectory(data_list)

set.seed(1234)
idx_train <- sample(1:L, N)
train_list <- data_list[idx_train]
test_list <- data_list[setdiff(1:L, idx_train)]
plot_trajectory(train_list)

train_mat <- create_features_matrix(train_list)
test_mat <- create_features_matrix(test_list)

############ 2) Set parameters for the 3 models ###############

########### 2.1) Set parameters for BB (Poisson and Neg-Bin)

# Set prior hyperparameters specific for poisson
a_l_poiss <- 50
b_l_poiss <- 0.1
print(paste0("E(lambda) = ", a_l_poiss/ b_l_poiss))
print(paste0("Var(lambda) = ", a_l_poiss/ (b_l_poiss^2)))

# Set prior hyperparameters specific for NB
q_star_nb <- 0.002
print(paste0("E(nstar) = ", 1/ q_star_nb))
print(paste0("Var(nstar) = ", (1-q_star_nb)/ (q_star_nb^2)) )
alpha_p_nb <- 0.000001
beta_p_nb <- 0.000001
print(paste0("E(p) = ", alpha_p_nb/(alpha_p_nb + beta_p_nb)))
print(paste0("Var(p) = ",(alpha_p_nb*beta_p_nb)/((alpha_p_nb + beta_p_nb)^2 *(alpha_p_nb+beta_p_nb+1)) ))

# Set other hyperparameters
a_alpha_bb <- 1
b_alpha_bb <- 0.1
print(paste0("E(alpha_bar) = ", a_alpha_bb/ b_alpha_bb))
print(paste0("Var(alpha_bar) = ", a_alpha_bb/ (b_alpha_bb^2)))
a_s_bb <- 1
b_s_bb <- 0.1
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


########## 2.2) Set parameters for the Gamma IBP

# Set prior hyperparameters
p_ibp <- 0.05
print(paste0("E(a) = ", 1/ p_ibp))
print(paste0("Var(a) = ", (1-p_ibp)/ (p_ibp^2)) )
r_ibp <- 1
t_ibp <- 0.1
print(paste0("E(b) = ", r_ibp/ t_ibp))
print(paste0("Var(b) = ", r_ibp/ (t_ibp^2)))
a_alpha_ibp <- 0.1
b_alpha_ibp <- 0.1
print(paste0("E(alpha) = ", a_alpha_ibp/ (a_alpha_ibp + b_alpha_ibp)))
print(paste0("Var(alpha) = ", a_alpha_ibp*b_alpha_ibp/ (a_alpha_ibp + b_alpha_ibp)^2 /(a_alpha_ibp+b_alpha_ibp+1)))
a_s_ibp <- 1
b_s_ibp <- 0.05
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


########## 2.3) Set MCMC parameters (common to all 3 models)

S_poiss <- S_negbin <- 4*10^4
S_ibp <- 4*10^4
n_burnin_poiss <- n_burnin_negbin <- 5*10^3
n_burnin_ibp <- 5*10^3
thin_poiss <- thin_negbin <- thin_ibp <- 5
seed <- 1234


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

############# 3.1) Processing output: BB + Poisson

# Mixing checks for (lambda, alpha, theta)
plot(lambda_chain_poiss, type="l") # mixing of lambda in prior number of features
# Running mean
plot(cumsum(lambda_chain_poiss)/(1:n_saved_iter_poiss), type="l")

plot(alpha_chain_poiss, type="l") # mixing of "alpha" 
# Running mean
plot(cumsum(alpha_chain_poiss)/(1:n_saved_iter_poiss), type="l")

plot(theta_chain_poiss, type="l") # mixing of "theta" in the beta prior on weights
# Running mean
plot(cumsum(theta_chain_poiss)/(1:n_saved_iter_poiss), type="l")


######## 3.2) Model-checking on Kn within sample: BB + Poisson

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


####### 4) Run BB + Negative-Binomial ################

# Set tau for MALA
tau_nb <- 0.0001

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


############# 4.1) Processing outputs: BB + Neg-Bin

# Mixing checks for (nstar, p, alpha, theta)
plot(nstar_chain_negbin, type="l") # mixing of nstar in prior number of features
# Running mean
plot(cumsum(nstar_chain_negbin)/(1:n_saved_iter_negbin), type="l")

plot(p_chain_negbin, type="l") # mixing of p in prior number of featuresplot(alpha_chain, type="l") # mixing of "alpha" 
# Running mean
plot(cumsum(p_chain_negbin)/(1:n_saved_iter_negbin), type="l")

plot(theta_chain_negbin, type="l") # mixing of "theta" in the beta prior on weights
# Running mean
plot(cumsum(theta_chain_negbin)/(1:n_saved_iter_negbin), type="l")

plot(alpha_chain_negbin, type="l") # mixing of "alpha" 
# Running mean
plot(cumsum(alpha_chain_negbin)/(1:n_saved_iter_negbin), type="l")


######## 4.2) Model-checking on Kn within sample: BB + Neg-Bin

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


############ 5) Run IBP + Gamma #################

# Set tau for MALA
tau_ibp <- 0.001

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

############ 5.1) Processing outputs: IBP + Gamma

# Mixing checks for (a, b, alpha, theta)
plot(a_chain_ibp, type="l") # mixing of a in prior number of features
# Running mean
plot(cumsum(a_chain_ibp)/(1:n_saved_iter_ibp), type="l")

plot(b_chain_ibp, type="l") # mixing of b in prior number of features
# Running mean
plot(cumsum(b_chain_ibp)/(1:n_saved_iter_ibp), type="l")

plot(theta_chain_ibp, type="l") # mixing of "theta" in the beta prior on weights
# Running mean
plot(cumsum(theta_chain_ibp)/(1:n_saved_iter_ibp), type="l")

plot(alpha_chain_ibp, type="l") # mixing of "alpha" 
# Running mean
plot(cumsum(alpha_chain_ibp)/(1:n_saved_iter_ibp), type="l")


######## 5.2) Model-checking on Kn within sample: IBP + Gamma

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


####### 6) Comparison between the Poisson, Neg-Bin and IBP + Gamma ##########

# Running mean alpha (compare only Poisson and Neg-Bin)
plot(cumsum(alpha_chain_poiss)/(1:n_saved_iter_poiss), type="l", col="red", ylab = "alpha")
lines(cumsum(alpha_chain_negbin)/(1:n_saved_iter_negbin), type="l", col="blue")

# Running mean theta (compare only Poisson and Neg-Bin)
plot(cumsum(theta_chain_poiss)/(1:n_saved_iter_poiss), type="l", col="red", ylab = "theta")
lines(cumsum(theta_chain_negbin)/(1:n_saved_iter_negbin), type="l", col="blue")

######### 6.1) Comparison: Model checking
gg_kn_poiss <- plot_Kn_median_and_sample(data_list = data_list, est_ci_poiss)
gg_kn_poiss <- gg_kn_poiss + ggtitle("Poiss")
gg_kn_negbin <- plot_Kn_median_and_sample(data_list = data_list, est_ci_negbin)
gg_kn_negbin <- gg_kn_negbin + ggtitle("Neg-Bin")
gg_kn_ibp <- plot_Kn_median_and_sample(data_list = data_list, est_ci_ibp)
gg_kn_ibp <- gg_kn_ibp + ggtitle("Gamma IBP")

ggarrange(gg_kn_poiss, gg_kn_negbin, gg_kn_ibp,
          ncol = 3, nrow = 1)

######## 6.2) Model checking with rarefaction
gg_kn_rar_poiss <- plot_Kn_median_and_rarefaction(train_list = train_list, est_ci_poiss, n_avg = 100)
gg_kn_rar_poiss <- gg_kn_rar_poiss + ggtitle("Poiss")
gg_kn_rar_negbin <- plot_Kn_median_and_rarefaction(train_list = train_list, est_ci_negbin, n_avg = 100)
gg_kn_rar_negbin <- gg_kn_rar_negbin + ggtitle("Neg-Bin")
gg_kn_rar_ibp <- plot_Kn_median_and_rarefaction(train_list = train_list, est_ci_ibp, n_avg = 100)
gg_kn_rar_ibp <- gg_kn_rar_ibp + ggtitle("Gamma IBP")

ggarrange(gg_kn_rar_poiss, gg_kn_rar_negbin, gg_kn_rar_ibp,
          ncol = 3, nrow = 1)

########### 6.3) Prediction on additional sample (plot for both training/test split and only training)
Kn = ncol(train_mat[,colSums(train_mat) > 0])

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


# Plot only training and prediction
gg_kmn_pred_poiss <- plot_Kmn_median_pred_and_rarefaction(train_list = train_list, 
                                                          ci = est_ci_pred_poiss, n_avg = 100)
gg_kmn_pred_poiss <- gg_kmn_pred_poiss + ggtitle("Poiss")
gg_kmn_pred_negbin <- plot_Kmn_median_pred_and_rarefaction(train_list = train_list, 
                                                           ci = est_ci_pred_negbin, n_avg = 100)
gg_kmn_pred_negbin <- gg_kmn_pred_negbin + ggtitle("Neg-Bin")
gg_kmn_pred_ibp <- plot_Kmn_median_pred_and_rarefaction(train_list = train_list, 
                                                        ci = est_ci_pred_ibp, n_avg = 100)
gg_kmn_pred_ibp <- gg_kmn_pred_ibp + ggtitle("Gamma IBP")

ggarrange(gg_kmn_pred_poiss, gg_kmn_pred_negbin, gg_kmn_pred_ibp,
          ncol = 3, nrow = 1)

# Plot training, prediction and test
gg_kmn_pred_test_poiss <- plot_Kmn_median_pred_and_test(train_list = train_list, test_list = test_list,
                                                        ci = est_ci_pred_poiss, n_avg = 100)
gg_kmn_pred_test_poiss <- gg_kmn_pred_test_poiss + ggtitle("Poiss")
gg_kmn_pred_test_negbin <- plot_Kmn_median_pred_and_test(train_list = train_list, test_list = test_list,
                                                         ci = est_ci_pred_negbin, n_avg = 100)
gg_kmn_pred_test_negbin <- gg_kmn_pred_test_negbin + ggtitle("Neg-Bin")
gg_kmn_pred_test_ibp <- plot_Kmn_median_pred_and_test(train_list = train_list, test_list = test_list,
                                                      ci = est_ci_pred_ibp, n_avg = 100)
gg_kmn_pred_test_ibp <- gg_kmn_pred_test_ibp + ggtitle("Gamma IBP")

ggarrange(gg_kmn_pred_test_poiss, gg_kmn_pred_test_negbin, gg_kmn_pred_test_ibp,
          ncol = 3, nrow = 1)




#############################################################
################# SETTING 2: ZIPF (asymptotic behaviour ) ########################
#############################################################

############ 1) Generate the data #####################

# set number of individuals
L <- 3000 # 1500
N <- 200 # dimension of the training (used as whole sample when we don't consider test/train)
M <- L - N
# set maximum number of features in the population
H <- 100

# value of xi
xi = 1.4

seed = 1234
# generate data from zipf, in form of binary matrix
data_mat <- rzipf(n = L, K = H, xi = xi, seed = seed)

# convert the binary matrix into list of features
data_list <- create_features_list(data_mat)
plot_trajectory(data_list)

set.seed(1234)
idx_train <- sample(1:L, N)
train_list <- data_list[idx_train]
test_list <- data_list[setdiff(1:L, idx_train)]
plot_trajectory(train_list)

train_mat <- create_features_matrix(train_list)
test_mat <- create_features_matrix(test_list)

############ 2) Set parameters for the 3 models ###############

########### 2.1) Set parameters for BB (Poisson and Neg-Bin)

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
a_alpha_bb <- 10
b_alpha_bb <- 1
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


########## 2.2) Set parameters for the Gamma IBP

# Set prior hyperparameters
p_ibp <- 0.05
print(paste0("E(a) = ", 1/ p_ibp))
print(paste0("Var(a) = ", (1-p_ibp)/ (p_ibp^2)) )
r_ibp <- 1
t_ibp <- 0.1
print(paste0("E(b) = ", r_ibp/ t_ibp))
print(paste0("Var(b) = ", r_ibp/ (t_ibp^2)))
a_alpha_ibp <- 0.5
b_alpha_ibp <- 0.5
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


########## 2.3) Set MCMC parameters (common to all 3 models)

S_poiss <- S_negbin <- 4*10^4
S_ibp <- 4*10^4
n_burnin_poiss <- n_burnin_negbin <- 5*10^3
n_burnin_ibp <- 5*10^3
thin_poiss <- thin_negbin <- thin_ibp <- 5
seed <- 1234


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

############# 3.1) Processing output: BB + Poisson

# Mixing checks for (lambda, alpha, theta)
plot(lambda_chain_poiss, type="l") # mixing of lambda in prior number of features
# Running mean
plot(cumsum(lambda_chain_poiss)/(1:n_saved_iter_poiss), type="l")

plot(alpha_chain_poiss, type="l") # mixing of "alpha" 
# Running mean
plot(cumsum(alpha_chain_poiss)/(1:n_saved_iter_poiss), type="l")

plot(theta_chain_poiss, type="l") # mixing of "theta" in the beta prior on weights
# Running mean
plot(cumsum(theta_chain_poiss)/(1:n_saved_iter_poiss), type="l")


######## 3.2) Model-checking on Kn within sample: BB + Poisson

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


####### 4) Run BB + Negative-Binomial ################

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


############# 4.1) Processing outputs: BB + Neg-Bin

# Mixing checks for (nstar, p, alpha, theta)
plot(nstar_chain_negbin, type="l") # mixing of nstar in prior number of features
# Running mean
plot(cumsum(nstar_chain_negbin)/(1:n_saved_iter_negbin), type="l")

plot(p_chain_negbin, type="l") # mixing of p in prior number of featuresplot(alpha_chain, type="l") # mixing of "alpha" 
# Running mean
plot(cumsum(p_chain_negbin)/(1:n_saved_iter_negbin), type="l")

plot(theta_chain_negbin, type="l") # mixing of "theta" in the beta prior on weights
# Running mean
plot(cumsum(theta_chain_negbin)/(1:n_saved_iter_negbin), type="l")

plot(alpha_chain_negbin, type="l") # mixing of "alpha" 
# Running mean
plot(cumsum(alpha_chain_negbin)/(1:n_saved_iter_negbin), type="l")


######## 4.2) Model-checking on Kn within sample: BB + Neg-Bin

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


############ 5) Run IBP + Gamma #################

# Set tau for MALA
tau_ibp <- 0.001

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

############ 5.1) Processing outputs: IBP + Gamma

# Mixing checks for (a, b, alpha, theta)
plot(a_chain_ibp, type="l") # mixing of a in prior number of features
# Running mean
plot(cumsum(a_chain_ibp)/(1:n_saved_iter_ibp), type="l")

plot(b_chain_ibp, type="l") # mixing of b in prior number of features
# Running mean
plot(cumsum(b_chain_ibp)/(1:n_saved_iter_ibp), type="l")

plot(theta_chain_ibp, type="l") # mixing of "theta" in the beta prior on weights
# Running mean
plot(cumsum(theta_chain_ibp)/(1:n_saved_iter_ibp), type="l")

plot(alpha_chain_ibp, type="l") # mixing of "alpha" 
# Running mean
plot(cumsum(alpha_chain_ibp)/(1:n_saved_iter_ibp), type="l")


######## 5.2) Model-checking on Kn within sample: IBP + Gamma

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


####### 6) Comparison between the Poisson, Neg-Bin and IBP + Gamma ##########

# Running mean alpha (compare only Poisson and Neg-Bin)
plot(cumsum(alpha_chain_poiss)/(1:n_saved_iter_poiss), type="l", col="red", ylab = "alpha")
lines(cumsum(alpha_chain_negbin)/(1:n_saved_iter_negbin), type="l", col="blue")

# Running mean theta (compare only Poisson and Neg-Bin)
plot(cumsum(theta_chain_poiss)/(1:n_saved_iter_poiss), type="l", col="red", ylab = "theta")
lines(cumsum(theta_chain_negbin)/(1:n_saved_iter_negbin), type="l", col="blue")

######### 6.1) Comparison: Model checking
gg_kn_poiss <- plot_Kn_median_and_sample(data_list = data_list, est_ci_poiss)
gg_kn_poiss <- gg_kn_poiss + ggtitle("Poiss")
gg_kn_negbin <- plot_Kn_median_and_sample(data_list = data_list, est_ci_negbin)
gg_kn_negbin <- gg_kn_negbin + ggtitle("Neg-Bin")
gg_kn_ibp <- plot_Kn_median_and_sample(data_list = data_list, est_ci_ibp)
gg_kn_ibp <- gg_kn_ibp + ggtitle("Gamma IBP")

ggarrange(gg_kn_poiss, gg_kn_negbin, gg_kn_ibp,
          ncol = 3, nrow = 1)

######## 6.2) Model checking with rarefaction
gg_kn_rar_poiss <- plot_Kn_median_and_rarefaction(train_list = train_list, est_ci_poiss, n_avg = 100)
gg_kn_rar_poiss <- gg_kn_rar_poiss + ggtitle("Poiss")
gg_kn_rar_negbin <- plot_Kn_median_and_rarefaction(train_list = train_list, est_ci_negbin, n_avg = 100)
gg_kn_rar_negbin <- gg_kn_rar_negbin + ggtitle("Neg-Bin")
gg_kn_rar_ibp <- plot_Kn_median_and_rarefaction(train_list = train_list, est_ci_ibp, n_avg = 100)
gg_kn_rar_ibp <- gg_kn_rar_ibp + ggtitle("Gamma IBP")

ggarrange(gg_kn_rar_poiss, gg_kn_rar_negbin, gg_kn_rar_ibp,
          ncol = 3, nrow = 1)

########### 6.3) Prediction on additional sample (plot for both training/test split and only training)
Kn = ncol(train_mat[,colSums(train_mat) > 0])

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


# Plot only training and prediction
gg_kmn_pred_poiss <- plot_Kmn_median_pred_and_rarefaction(train_list = train_list, 
                                                          ci = est_ci_pred_poiss, n_avg = 100)
gg_kmn_pred_poiss <- gg_kmn_pred_poiss + ggtitle("Poiss")
gg_kmn_pred_negbin <- plot_Kmn_median_pred_and_rarefaction(train_list = train_list, 
                                                           ci = est_ci_pred_negbin, n_avg = 100)
gg_kmn_pred_negbin <- gg_kmn_pred_negbin + ggtitle("Neg-Bin")
gg_kmn_pred_ibp <- plot_Kmn_median_pred_and_rarefaction(train_list = train_list, 
                                                        ci = est_ci_pred_ibp, n_avg = 100)
gg_kmn_pred_ibp <- gg_kmn_pred_ibp + ggtitle("Gamma IBP")

ggarrange(gg_kmn_pred_poiss, gg_kmn_pred_negbin, gg_kmn_pred_ibp,
          ncol = 3, nrow = 1)

# Plot training, prediction and test
gg_kmn_pred_test_poiss <- plot_Kmn_median_pred_and_test(train_list = train_list, test_list = test_list,
                                                        ci = est_ci_pred_poiss, n_avg = 100)
gg_kmn_pred_test_poiss <- gg_kmn_pred_test_poiss + ggtitle("Poiss")
gg_kmn_pred_test_negbin <- plot_Kmn_median_pred_and_test(train_list = train_list, test_list = test_list,
                                                         ci = est_ci_pred_negbin, n_avg = 100)
gg_kmn_pred_test_negbin <- gg_kmn_pred_test_negbin + ggtitle("Neg-Bin")

ggarrange(gg_kmn_pred_test_poiss, gg_kmn_pred_test_negbin, 
          ncol = 2, nrow = 1)


gg_kmn_pred_test_ibp <- plot_Kmn_median_pred_and_test(train_list = train_list, test_list = test_list,
                                                      ci = est_ci_pred_ibp, n_avg = 100)
gg_kmn_pred_test_ibp <- gg_kmn_pred_test_ibp + ggtitle("Gamma IBP")

ggarrange(gg_kmn_pred_test_poiss, gg_kmn_pred_test_negbin, gg_kmn_pred_test_ibp,
          ncol = 3, nrow = 1)

# Limit distributions for Ntilde
Kn = ncol(train_mat[,colSums(train_mat) > 0])

ntilde_chain_poiss <- generate_Ntilde_chain_poiss(lambda_chain_poiss, alpha_chain_poiss, 
                                                  theta_chain_poiss, n = N, Kn)
ntilde_chain_poiss <- as.data.frame(ntilde_chain_poiss)

ntilde_chain_negbin <- generate_Ntilde_chain_negbin(nstar_chain_negbin, p_chain_negbin,
                                                    alpha_chain_negbin,
                                                    theta_chain_negbin, n = N, Kn)
ntilde_chain_negbin <- as.data.frame(ntilde_chain_negbin)

ggplot(ntilde_chain_poiss, aes(x=ntilde_chain_poiss)) + geom_density()
ggplot(ntilde_chain_negbin, aes(x=ntilde_chain_negbin)) + geom_density()
