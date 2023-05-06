rm(list=ls())
seed = 1234

# set true parameters of the poiss BB process
alpha <- -1
theta <- 10
nstar <- 10
p <- 0.05
# compute the corresponding mean and variance of the negbin
mu <- nstar*(1-p)/p
sigma2 <- nstar*(1-p)/(p**2)

# total number of samples
L = 5000

# generate data in list form
data_list_full <- buffet_negbin_BB(alpha = alpha, theta = theta, 
                                   n = L, nstar = nstar, p = p)
data_list <- data_list_full$features


data_mat <- create_features_matrix(data_list)

ncol(data_mat)

# Set prior hyperparameters
a_star <- 1
b_star <- 0.1
print(paste0("E(nstar) = ", a_star/ b_star))
print(paste0("Var(nstar) = ", a_star/ (b_star^2)))
alpha_p <- 0.000001
beta_p <- 0.000001
print(paste0("E(p) = ", alpha_p/(alpha_p + beta_p)))
print(paste0("Var(p) = ",(alpha_p*beta_p)/((alpha_p + beta_p)^2 *(alpha_p+beta_p+1)) ))
a_alpha <- 10
b_alpha <- 1
print(paste0("E(alpha_bar) = ", a_alpha/ b_alpha))
print(paste0("Var(alpha_bar) = ", a_alpha/ (b_alpha^2)))
a_s <- 2
b_s <- 0.1
print(paste0("E(s) = ", a_s/ b_s))
print(paste0("Var(s) = ", a_s/ (b_s^2)))

print(paste0("E(theta) = ", a_s/ b_s + a_alpha/ b_alpha))
print(paste0("Var(theta) = ", a_s/ (b_s^2) + a_alpha/ (b_alpha^2)))


# Set initial values for the parameters
nstar_0 <- 10
p_0 <- 0.8
alpha_bar_0 <- 1
s_0 <- 1

# Set tau for MALA
# tau <- 0.001 # when nstar is fixed
tau <- 0.00001

# Set MCMC parameters
S <- 5*10^5
n_burnin <- 0
thin <- 10
seed <- 1234

output <- gibbs_sampler_negbin(Z = data_mat,
                               nstar_0, p_0 = 0.01,  s_0 = 9 , alpha_bar_0 = 1,
                               a_star, b_star, alpha_p, beta_p, a_alpha, b_alpha, a_s, b_s,
                               tau, fixed = c(F, T,T,T),
                               S, n_burnin, thin, seed)

n_saved_iter <- length(output$nstar_vec)
nstar_chain <- output$nstar_vec
p_chain <- output$p_vec
s_chain <- output$s_vec
alpha_bar_chain <- output$alpha_bar_vec
alpha_chain <- - alpha_bar_chain
theta_chain <- s_chain + alpha_bar_chain

################################################################
############# Processing outputs ##############################
################################################################

# Mixing checks for (nstar, p, s, alpha_bar)
plot(nstar_chain, type="l") # mixing of nstar in prior number of features
plot(p_chain, type="l") # mixing of p in prior number of features
plot(s_chain, type="l") # mixing of "s" 
plot(alpha_bar_chain, type="l") # mixing of "alpha_bar" 

# Mixing checks for (nstar, p, alpha, theta)
plot(nstar_chain[4000:n_saved_iter], type="l") # mixing of nstar in prior number of features
plot(p_chain, type="l") # mixing of p in prior number of featuresplot(alpha_chain, type="l") # mixing of "alpha" 
plot(theta_chain, type="l") # mixing of "theta" in the beta prior on weights
plot(alpha_chain, type="l") # mixing of "alpha" 

alpha 
theta 
nstar
p 

