rm(list=ls())
seed = 1234

# set true parameters of the poiss BB process
alpha_true <- 0.8
theta_true <- 50
a_true <- 20
b_true <- 10
s_true <- theta_true + alpha_true

# total number of samples
L = 300

# generate data in list form
data_list_full <- buffet_gamma_IBP(alpha = alpha_true, theta = theta_true, 
                                   n = L, a = a_true, b = b_true)
data_list <- data_list_full$features


data_mat <- create_features_matrix(data_list)

ncol(data_mat)

plot_trajectory(data_list)


##########################################################################

# Set prior hyperparameters
p <- 0.05
print(paste0("E(a) = ", 1/ p))
print(paste0("Var(a) = ", (1-p)/ (p^2)) )
r <- 0.1
t <- 0.01
print(paste0("E(b) = ", r/ t))
print(paste0("Var(b) = ", r/ (t^2)))
a_alpha <- 0.1
b_alpha <- 0.1
print(paste0("E(alpha) = ", a_alpha/ (a_alpha + b_alpha)))
print(paste0("Var(alpha) = ", a_alpha*b_alpha/ (a_alpha + b_alpha)^2 /(a_alpha+b_alpha+1)))
a_s <- 2
b_s <- 0.05
print(paste0("E(s) = ", a_s/ b_s))
print(paste0("Var(s) = ", a_s/ (b_s^2)))

print(paste0("E(theta) = ", a_s/ b_s -  a_alpha/ (a_alpha + b_alpha)))
print(paste0("Var(theta) = ", a_s/ (b_s^2) + a_alpha*b_alpha/ (a_alpha + b_alpha)^2 /(a_alpha+b_alpha+1)))


# Set initial values for the parameters
a_0 <- 5
b_0 <- 1
alpha_0 <- 0.5
s_0 <- 15

# Set tau for MALA
tau <- 0.001

# Set MCMC parameters
S <- 4*10^4
n_burnin <- 5*10^3
thin <- 5
seed <- 1234

output <- gibbs_sampler_gamma_ibp(Z = data_mat,
                                  a_0 = a_true, b_0 = b_true,  s_0 = s_true, alpha_0 = alpha_true,
                                  p, r, t, a_alpha, b_alpha, a_s, b_s,
                                  tau, fixed = c(F,F,F,F),
                                  S, n_burnin, thin, seed)

n_saved_iter <- length(output$a_vec)
a_chain <- output$a_vec
b_chain <- output$b_vec
s_chain <- output$s_vec
alpha_chain <- output$alpha_vec
theta_chain <- s_chain - alpha_chain

################################################################
############# Processing outputs ##############################
################################################################

# Mixing checks for (nstar, p, s, alpha_bar)
plot(a_chain, type="l") # mixing of a in prior number of features
plot(b_chain, type="l") # mixing of b in prior number of features
plot(s_chain, type="l") # mixing of "s" 
plot(alpha_chain, type="l") # mixing of "alpha" 

# Mixing checks for (a, b, alpha, theta)
plot(a_chain, type="l") # mixing of a in prior number of features
# Running mean
plot(cumsum(a_chain)/(1:n_saved_iter), type="l")

plot(b_chain, type="l") # mixing of b in prior number of features
# Running mean
plot(cumsum(b_chain)/(1:n_saved_iter), type="l")

plot(theta_chain, type="l") # mixing of "theta" in the beta prior on weights
# Running mean
plot(cumsum(theta_chain)/(1:n_saved_iter), type="l")

plot(alpha_chain, type="l") # mixing of "alpha" 
# Running mean
plot(cumsum(alpha_chain)/(1:n_saved_iter), type="l")

##############################################################
######## Model-checking on Kn within sample ##################
##############################################################
M <- nrow(data_mat)
kmn_chain <- generate_Kmn_chain_gamma_ibp(a_chain, b_chain, alpha_chain, 
                                       theta_chain, M, n=0)

est_ci_gamma_ibp <- matrix(NA, nrow = M, ncol = 3)
# first column = lower bound
# second columns = medians
# third columns = upper bound
for (m in 1:M){
  est_ci_gamma_ibp[m,] <- quantile(kmn_chain[m,], probs = c(0.025,0.5,0.975))
}
est_ci_gamma_ibp <- list("medians" = est_ci_gamma_ibp[,2],
                      "lbs" = est_ci_gamma_ibp[,1],
                      "ubs" = est_ci_gamma_ibp[,3])

plot_Kn_median_and_sample(data_list = data_list, est_ci_gamma_ibp)

