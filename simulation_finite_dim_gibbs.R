##############################################################
###### Simulation in FINITE-DIMENSIONAL CASE #################
##############################################################
rm(list=ls())
library(ggpubr)

#############################################################
################# SETTING 1: UNIFORM ########################
#############################################################

# set number of individuals
L <- 1000
N <- 600 # dimension of the training (used as whole sample when we don't consider test/train)
M <- L - N
# set maximum number of features in the population
H <- 300

# set probability of presence of each feature in each individual
prob_presence <- 0.005

set.seed(1234)
data_mat <- matrix(rbinom(L*H, size = 1, prob = prob_presence), nrow = L, ncol = H)

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

############ BB + Poisson and BB + Negative-Binomial #####################

# Set prior hyperparameters specific for poisson
a_l <- 50
b_l <- 0.1
print(paste0("E(lambda) = ", a_l/ b_l))
print(paste0("Var(lambda) = ", a_l/ (b_l^2)))
# Set prior hyperparameters specific for NB
q_star <- 0.002
print(paste0("E(nstar) = ", 1/ q_star))
print(paste0("Var(nstar) = ", (1-q_star)/ (q_star^2)) )
alpha_p <- 0.000001
beta_p <- 0.000001
print(paste0("E(p) = ", alpha_p/(alpha_p + beta_p)))
print(paste0("Var(p) = ",(alpha_p*beta_p)/((alpha_p + beta_p)^2 *(alpha_p+beta_p+1)) ))
# Set other hyperparameters
a_alpha <- 1
b_alpha <- 0.1
print(paste0("E(alpha_bar) = ", a_alpha/ b_alpha))
print(paste0("Var(alpha_bar) = ", a_alpha/ (b_alpha^2)))
a_s <- 1
b_s <- 0.1
print(paste0("E(s) = ", a_s/ b_s))
print(paste0("Var(s) = ", a_s/ (b_s^2)))

# Set initial values for poisson 
lambda_0 <- 100
# Set initial values for NB
nstar_0 <- 100
p_0 <- 0.2
# Set initial values for other parameters
alpha_bar_0 <- 1
s_0 <- 1


# Set MCMC parameters
S_poiss <- S_negbin <- 5*10^4
n_burnin_poiss <- n_burnin_negbin <- 5*10^3
thin_poiss <- thin_negbin <- 5
seed <- 1234

########### Run BB + Poisson 

# Set tau for MALA
tau <- 0.0005

output_poiss <- gibbs_sampler_poiss(Z = train_mat,
                              lambda_0, alpha_bar_0, s_0,
                              a_l, b_l, a_alpha, b_alpha, a_s, b_s,
                              tau,
                              S_poiss, n_burnin_poiss, thin_poiss, seed)

n_saved_iter_poiss <- length(output_poiss$lambda_vec)
lambda_chain_poiss <- output_poiss$lambda_vec
s_chain_poiss <- output_poiss$s_vec
alpha_bar_chain_poiss <- output_poiss$alpha_bar_vec
alpha_chain_poiss <- - alpha_bar_chain_poiss
theta_chain_poiss <- s_chain_poiss + alpha_bar_chain_poiss

################################################################
############# Processing outputs ##############################
################################################################

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

##############################################################
######## Model-checking on Kn within sample ##################
##############################################################
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


####### Run BB + Negative-Binomial

# Set tau for MALA
tau <- 0.0002

output_negbin <- gibbs_sampler_negbin_geometric(Z = train_mat,
                                         nstar_0 = 100, p_0 = 0.2,  s_0 = 20, alpha_bar_0 = 5,
                                         q_star, alpha_p, beta_p, a_alpha, b_alpha, a_s, b_s,
                                         tau, fixed = c(F,F,F,F),
                                         S_negbin, n_burnin_negbin, thin_negbin, seed)

n_saved_iter_negbin <- length(output_negbin$nstar_vec)
nstar_chain_negbin <- output_negbin$nstar_vec
p_chain_negbin <- output_negbin$p_vec
s_chain_negbin <- output_negbin$s_vec
alpha_bar_chain_negbin <- output_negbin$alpha_bar_vec
alpha_chain_negbin <- - alpha_bar_chain_negbin
theta_chain_negbin <- s_chain_negbin + alpha_bar_chain_negbin

################################################################
############# Processing outputs ##############################
################################################################

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

##############################################################
######## Model-checking on Kn within sample ##################
##############################################################
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


######## Comparison between the Poisson and NegBin
# Running mean alpha
plot(cumsum(alpha_chain_poiss)/(1:n_saved_iter_poiss), type="l", col="red", ylab = "alpha")
lines(cumsum(alpha_chain_negbin)/(1:n_saved_iter_negbin), type="l", col="blue")

# Running mean theta
plot(cumsum(theta_chain_poiss)/(1:n_saved_iter_poiss), type="l", col="red", ylab = "theta")
lines(cumsum(theta_chain_negbin)/(1:n_saved_iter_negbin), type="l", col="blue")

# Model checking
gg_kn_poiss <- plot_Kn_median_and_sample(data_list = data_list, est_ci_poiss)
gg_kn_poiss <- gg_kn_poiss + ggtitle("Poiss")
gg_kn_negbin <- plot_Kn_median_and_sample(data_list = data_list, est_ci_negbin)
gg_kn_negbin <- gg_kn_negbin + ggtitle("Neg-Bin")

ggarrange(gg_kn_poiss, gg_kn_negbin,
          ncol = 2, nrow = 1)

# Model checking with rarefaction
gg_kn_rar_poiss <- plot_Kn_median_and_rarefaction(train_list = train_list, est_ci_poiss, n_avg = 100)
gg_kn_rar_poiss <- gg_kn_rar_poiss + ggtitle("Poiss")
gg_kn_rar_negbin <- plot_Kn_median_and_rarefaction(train_list = train_list, est_ci_negbin, n_avg = 100)
gg_kn_rar_negbin <- gg_kn_rar_negbin + ggtitle("Neg-Bin")

ggarrange(gg_kn_rar_poiss, gg_kn_rar_negbin,
          ncol = 2, nrow = 1)

########### Prediction on additional sample (plot for both training/test split and only training)
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
                                             M = M, n = N)

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

# Plot only training and prediction
gg_kmn_pred_poiss <- plot_Kmn_median_pred_and_rarefaction(train_list = train_list, 
                                                          ci = est_ci_pred_poiss, n_avg = 100)
gg_kmn_pred_poiss <- gg_kmn_pred_poiss + ggtitle("Poiss")
gg_kmn_pred_negbin <- plot_Kmn_median_pred_and_rarefaction(train_list = train_list, 
                                                          ci = est_ci_pred_negbin, n_avg = 100)
gg_kmn_pred_negbin <- gg_kmn_pred_negbin + ggtitle("Neg-Bin")

ggarrange(gg_kmn_pred_poiss, gg_kmn_pred_negbin,
          ncol = 2, nrow = 1)

# Plot training, prediction and test
gg_kmn_pred_test_poiss <- plot_Kmn_median_pred_and_test(train_list = train_list, test_list = test_list,
                                                          ci = est_ci_pred_poiss, n_avg = 100)
gg_kmn_pred_test_poiss <- gg_kmn_pred_test_poiss + ggtitle("Poiss")
gg_kmn_pred_test_negbin <- plot_Kmn_median_pred_and_test(train_list = train_list, test_list = test_list,
                                                           ci = est_ci_pred_negbin, n_avg = 100)
gg_kmn_pred_test_negbin <- gg_kmn_pred_test_negbin + ggtitle("Neg-Bin")

ggarrange(gg_kmn_pred_test_poiss, gg_kmn_pred_test_negbin,
          ncol = 2, nrow = 1)


###### Plot limiting distribution of Ntilde
# Poisson 
Ntilde_chain_poiss <- generate_Ntilde_chain_poiss(lambda_chain_poiss, alpha_chain_poiss,
                                            theta_chain_poiss, n = N)

# Negative Binomial
Ntilde_chain_negbin <- generate_Ntilde_chain_negbin(nstar_chain_negbin, p_chain_negbin,
                                              alpha_chain_negbin, theta_chain_negbin,
                                              n = N)
#############################################################
################# SETTING 2: ZIPF ###########################
#############################################################
