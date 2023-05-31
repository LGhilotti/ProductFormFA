####################################################################
## Simulation study using Poiss-BB model - known hyperparameters ###
####################################################################
library(ggpubr)

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
  
  # plot the E[KMN|Z_N] together with 90% credible intervals, for given N and train set
  CI_given_sample <- CI_Kmn_poiss_BB(alpha, theta, M, N, lambda, 0.9)
  ggf <- plot_Kmn_given_sample(N, train_list, CI_given_sample)
  
  gg_list[[j]] <- ggf
  
}

# display plots on the same page
ggarrange(gg_list[[1]], gg_list[[2]], gg_list[[3]],
          ncol = 3, nrow = 1)

####################################################################
## Simulation study using Poiss-BB model - unknown hyperparameters -
## EFPF-maximization 
####################################################################

seed = 123

# set true parameters of the poiss BB process
alpha <- -1
theta <- 10
lambda <- 1000

# total number of samples
L = 1000

# generate data in list form
data_list_full <- buffet_poiss_BB(alpha = alpha, theta = theta, 
                             n = L, lambda = lambda)
data_list <- data_list_full$features

# set the dimension of the training set and the training set itself
Ns <- c(100, 200, 400)

# list of ggplots
gg_list <- vector("list", length = length(Ns) )

for (j in 1:length(Ns)){
  N <- Ns[j]

  train_list <- data_list[1:N]
  
  # set the dimension of the test set and the test set itself
  M <- L - N
  test_list <- data_list[(N+1):L]
  
  # EB: fix the model hyperparameters based on EFPF-maximization
  train_counts <- tabulate(unlist(train_list))[tabulate(unlist(train_list)) != 0]
  eb_efpf_poiss_BB <- EB_EFPF_poiss_BB(n = N, 
                                       counts = train_counts, pars_0 = c(-2, 5, 10))
  
  eb_efpf_poiss_BB
  alpha_est <- eb_efpf_poiss_BB[1]
  theta_est <- eb_efpf_poiss_BB[2]
  lambda_est <- eb_efpf_poiss_BB[3]
  
  # plot the E[KMN|Z_N] together with 90% credible intervals, for given N and train set
  CI_given_sample <- CI_Kmn_poiss_BB(alpha_est, theta_est, M, N, lambda_est, 0.9)
  ggf <- plot_Kmn_given_sample_with_observed(N, data_list, CI_given_sample)
  
  gg_list[[j]] <- ggf
  
}

# display plots on the same page
ggarrange(gg_list[[1]], gg_list[[2]], gg_list[[3]],
          ncol = 3, nrow = 1)


####################################################################
## Simulation study using Poiss-BB model - unknown hyperparameters -
## Method of moments 
####################################################################

seed = 123

# set true parameters of the poiss BB process
alpha <- -1
theta <- 50
lambda <- 1000

# total number of samples
L = 1000

# generate data in list form
data_list_full <- buffet_poiss_BB(alpha = alpha, theta = theta, 
                                  n = L, lambda = lambda)
data_list <- data_list_full$features

# set the dimension of the training set and the training set itself
Ns <- c(100, 200, 400)

# list of ggplots
gg_list <- vector("list", length = length(Ns) )

for (j in 1:length(Ns)){
  N <- Ns[j]
  
  train_list <- data_list[1:N]
  
  # set the dimension of the test set and the test set itself
  M <- L - N
  test_list <- data_list[(N+1):L]
  
  # EB: fix the model hyperparameters based on Method of moments
  train_counts <- tabulate(unlist(train_list))[tabulate(unlist(train_list)) != 0]
  eb_mm_poiss_BB <- EB_MM_poiss_BB(n = N, 
                                   ntrain = round(N*2/3),
                                   val_rep = 5, data_list = data_list,
                                   pars_0 = c(-1,10,10), seed_id = 1234)
  
  alpha_est <- eb_mm_poiss_BB[1]
  theta_est <- eb_mm_poiss_BB[2]
  lambda_est <- eb_mm_poiss_BB[3]
  
  # plot the E[KMN|Z_N] together with 90% credible intervals, for given N and train set
  CI_given_sample <- CI_Kmn_poiss_BB(alpha_est, theta_est, M, N, lambda_est, 0.9)
  ggf <- plot_Kmn_given_sample_with_observed(N, data_list, CI_given_sample)
  
  gg_list[[j]] <- ggf
  
}

# display plots on the same page
ggarrange(gg_list[[1]], gg_list[[2]], gg_list[[3]],
          ncol = 3, nrow = 1)


