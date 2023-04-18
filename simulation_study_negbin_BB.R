####################################################################
## Simulation study using NegBin-BB model - known hyperparameters ###
####################################################################
library(ggpubr)

seed = 1234

# set true parameters of the poiss BB process
alpha <- -1
theta <- 10
nstar <- 1000
p <- 0.8
# compute the corresponding mean and variance of the negbin
mu <- nstar*(1-p)/p
sigma2 <- nstar*(1-p)/(p**2)

# total number of samples
L = 200

# generate data in list form
data_list_full <- buffet_negbin_BB(alpha = alpha, theta = theta, 
                                  n = L, nstar = nstar, p = p)
data_list <- data_list_full$features

# set the dimension of the training set and the training set itself
Ns <- c(10, 50, 100)

# list of ggplots
gg_list <- vector("list", length = length(Ns) )

for (j in 1:length(Ns)){
  N <- Ns[j]
  
  train_list <- data_list[1:N]
  
  # set the dimension of the test set and the test set itself
  M <- L - N
  test_list <- data_list[(N+1):L]
  
  # compute number of observed feature in the training set Kn
  Kn <- length(unique(unlist(train_list)))
  
  # plot the E[KMN|Z_N] together with 90% credible intervals, for given N and train set
  CI_given_sample <- CI_Kmn_negbin_BB(alpha, theta, M, N, Kn, nstar, p, 0.9)
  ggf <- plot_Kmn_given_sample(N, train_list, CI_given_sample)
  
  gg_list[[j]] <- ggf
  
}

# display plots on the same page
ggarrange(gg_list[[1]], gg_list[[2]], gg_list[[3]],
          ncol = 3, nrow = 1)

####################################################################
## Simulation study using NegBin-BB model - unknown hyperparameters -
## EFPF-maximization
####################################################################

seed = 12345

# set true parameters of the poiss BB process
alpha <- -1
theta <- 15
nstar <- 1000
p <- 0.8
# compute the corresponding mean and variance of the negbin
mu <- nstar*(1-p)/p
sigma2 <- nstar*(1-p)/(p**2)

# total number of samples
L = 1000

# generate data in list form
data_list_full <- buffet_negbin_BB(alpha = alpha, theta = theta, 
                                  n = L, nstar = nstar, p = p)
data_list <- data_list_full$features

# set the dimension of the training set and the training set itself
Ns <- c(50, 100, 150)

# list of ggplots
gg_list <- vector("list", length = length(Ns) )

for (j in 1:length(Ns)){
  N <- Ns[j]
  
  train_list <- data_list[1:N]
  
  # set the dimension of the test set and the test set itself
  M <- L - N
  test_list <- data_list[(N+1):L]
  
  # compute number of observed feature in the training set Kn
  Kn <- length(unique(unlist(train_list)))
  
  # EB: fix the model hyperparameters based on EFPF-maximization
  train_counts <- tabulate(unlist(train_list))[tabulate(unlist(train_list)) != 0]
  eb_efpf_negbin_BB <- EB_EFPF_fixed_negbin_mv_BB(n = N, 
                                       counts = train_counts, pars_0 = c(-2, 5, 100, 150))
  
  eb_efpf_negbin_BB
  alpha_est <- eb_efpf_negbin_BB[1]
  theta_est <- eb_efpf_negbin_BB[2]
  mu_est <- eb_efpf_negbin_BB[3]
  sigma2_est <- eb_efpf_negbin_BB[4]
  
  nstar_est <- (mu_est**2)/(sigma2_est - mu_est)
  p_est <- mu_est / sigma2_est
  
  # plot the E[KMN|Z_N] together with 90% credible intervals, for given N and train set
  CI_given_sample <- CI_Kmn_negbin_BB(alpha_est, theta_est, M, N, Kn, 
                                      nstar_est, p_est, 0.9)
  ggf <- plot_Kmn_given_sample_with_observed(N, data_list, CI_given_sample)
  
  gg_list[[j]] <- ggf
  
}

# display plots on the same page
ggarrange(gg_list[[1]], gg_list[[2]], gg_list[[3]],
          ncol = 3, nrow = 1)

####################################################################
## Simulation study using NegBin-BB model - unknown hyperparameters -
## Method of moments
####################################################################

seed = 1234

# set true parameters of the poiss BB process
alpha <- -1
theta <- 50
nstar <- 10000
p <- 0.5
# compute the corresponding mean and variance of the negbin
mu <- nstar*(1-p)/p
sigma2 <- nstar*(1-p)/(p**2)

# total number of samples
L = 1000

# generate data in list form
data_list_full <- buffet_negbin_BB(alpha = alpha, theta = theta, 
                                   n = L, nstar = nstar, p = p)
data_list <- data_list_full$features

# set the dimension of the training set and the training set itself
Ns <- c(50, 100, 150)

# list of ggplots
gg_list <- vector("list", length = length(Ns) )

for (j in 1:length(Ns)){
  N <- Ns[j]
  
  train_list <- data_list[1:N]
  
  # set the dimension of the test set and the test set itself
  M <- L - N
  test_list <- data_list[(N+1):L]
  
  # compute number of observed feature in the training set Kn
  Kn <- length(unique(unlist(train_list)))
  
  # EB: fix the model hyperparameters based on EFPF-maximization
  train_counts <- tabulate(unlist(train_list))[tabulate(unlist(train_list)) != 0]
  eb_mm_negbin_mv_BB <- EB_MM_negbin_mv_BB(n = N, 
                                        ntrain = round(N*2/3),
                                        val_rep = 5, data_list = data_list,
                                        pars_0 = c(-2,5,10,50), seed_id = 123)
  
  alpha_est <- eb_mm_negbin_mv_BB[1]
  theta_est <- eb_mm_negbin_mv_BB[2]
  mu_est <- eb_mm_negbin_mv_BB[3]
  sigma2_est <- eb_mm_negbin_mv_BB[4]
  
  nstar_est <- (mu_est**2)/(sigma2_est - mu_est)
  p_est <- mu_est / sigma2_est
  
  # plot the E[KMN|Z_N] together with 90% credible intervals, for given N and train set
  CI_given_sample <- CI_Kmn_negbin_BB(alpha_est, theta_est, M, N, Kn, 
                                      nstar_est, p_est, 0.9)
  ggf <- plot_Kmn_given_sample_with_observed(N, data_list, CI_given_sample)
  
  gg_list[[j]] <- ggf
  
}

# display plots on the same page
ggarrange(gg_list[[1]], gg_list[[2]], gg_list[[3]],
          ncol = 3, nrow = 1)

