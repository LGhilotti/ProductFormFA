####################################################################
#### Simulation study using Poiss-BB model #########################
####################################################################

seed <- 1234

# set true parameters of the poiss BB process
alpha <- -5
theta <- 15
lambda <- 1000

# total number of samples
L <- 50

# generate data in list form
data_list_full <- buffet_poiss_BB(
  alpha = alpha, theta = theta,
  n = L, lambda = lambda
)
data_list <- data_list_full$features

# set the dimension of the training set and the training set itself
N <- 10
train_list <- data_list[1:N]

# set the dimension of the test set and the test set itself
M <- L - N
test_list <- data_list[(N + 1):L]

# EB: fix the model hyperparameters based on EFPF-maximization
train_counts <- tabulate(unlist(train_list))
eb_efpf_poiss_BB <- EB_EFPF_poiss_BB(
  n = N,
  counts = train_counts, pars_0 = c(-2, 5, 10)
)

eb_efpf_poiss_BB
alpha_est <- eb_efpf_poiss_BB[1]
theta_est <- eb_efpf_poiss_BB[2]
lambda_est <- eb_efpf_poiss_BB[3]

# estimated number of new features
est_new_features <- mean_kmn_poiss_BB(
  alpha = alpha_est, theta = theta_est,
  m = M, n = N, lambda = lambda_est
)

# compute the percentage accuracy of the estimate with respect to the test set
perc_acc_list <- perc_accuracy(train_list, test_list, est_new_features)

# plot the E[KMN|Z_N] together with 90% credible intervals, for given N and train set
CI_given_sample <- CI_Kmn_poiss_BB(alpha_est, theta_est, M, N, lambda_est, 0.9)
plot_Kmn_given_sample_with_observed(N, data_list, CI_given_sample)
