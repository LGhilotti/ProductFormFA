rm(list=ls())
seed = 1234

# set true parameters of the poiss BB process
alpha <- -1
theta <- 10
lambda <- 1000

# total number of samples
L = 500

# generate data in list form
data_list_full <- buffet_poiss_BB(alpha = alpha, theta = theta, 
                                  n = L, lambda = lambda)
data_list <- data_list_full$features

data_mat <- create_features_matrix(data_list)


# Set prior hyperparameters
a_l <- 10
b_l <- 0.1
print(paste0("E(lambda) = ", a_l/ b_l))
print(paste0("Var(lambda) = ", a_l/ (b_l^2)))
a_alpha <- 10
b_alpha <- 1
print(paste0("E(alpha_bar) = ", a_alpha/ b_alpha))
print(paste0("Var(alpha_bar) = ", a_alpha/ (b_alpha^2)))
a_s <- 0.1
b_s <- 0.1
print(paste0("E(s) = ", a_s/ b_s))
print(paste0("Var(s) = ", a_s/ (b_s^2)))

# Set initial values for the parameters
lambda_0 <- 100
alpha_bar_0 <- 1
s_0 <- 1

# Set tau for MALA
tau <- 0.001

# Set MCMC parameters
S <- 10^4
n_burnin <- 10^3
thin <- 5
seed <- 1234

output <- gibbs_sampler_poiss(Z = data_mat,
                              lambda_0, alpha_bar_0, s_0,
                              a_l, b_l, a_alpha, b_alpha, a_s, b_s,
                              tau,
                              S, n_burnin, thin, seed)

n_saved_iter <- length(output$lambda_vec)
lambda_chain <- output$lambda_vec
s_chain <- output$s_vec
alpha_bar_chain <- output$alpha_bar_vec
alpha_chain <- - alpha_bar_chain
theta_chain <- s_chain + alpha_bar_chain

################################################################
############# Processing outputs ##############################
################################################################

# Mixing checks for (lambda, s, alpha_bar)
plot(lambda_chain, type="l") # mixing of lambda in prior number of features
plot(s_chain, type="l") # mixing of "s" 
plot(alpha_bar_chain, type="l") # mixing of "alpha_bar" 

# Mixing checks for (lambda, alpha, theta)
plot(lambda_chain, type="l") # mixing of lambda in prior number of features
plot(alpha_chain, type="l") # mixing of "alpha" 
plot(theta_chain, type="l") # mixing of "theta" in the beta prior on weights
