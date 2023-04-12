####################################################################
#### Simulation study using Zipf-distributed data ##################
####################################################################

seed = 1234

# total number of samples
L = 2000
# value of xi
xi = 1

# generate data from zipf, in form of binary matrix
data_mat <- rzipf(n = L, K = 10^5, xi = xi, seed = seed)

# convert the binary matrix into list of features
data_list <- create_features_list(data_mat)

# set the dimension of the training set and the training set itself
N <- 100
train_list <- data_list[1:N]
train_mat <- data_mat[1:N, ]

# set the dimension of the test set and the test set itself
M <- L - N
test_list <- data_list[(N+1):L]
test_mat <- data_mat[(N+1):L, ]

# EB: fix the model hyperparameters based on EFPF-maximization
train_counts <- colSums(train_mat)[colSums(train_mat)!=0]
eb_efpf_poiss_BB <- EB_EFPF_poiss_BB(n = N, 
                                     counts = train_counts, pars_0 = c(-5, 10, 100))

eb_efpf_poiss_BB
alpha_est <- eb_efpf_poiss_BB[1]
theta_est <- eb_efpf_poiss_BB[2]
lambda_est <- eb_efpf_poiss_BB[3]

# estimated number of new features
est_new_features <- mean_kmn_poiss_BB(alpha = alpha_est, theta = theta_est, 
                                      m = M, n = N, lambda = lambda_est)

# compute the percentage accuracy of the estimate with respect to the test set
perc_acc_list <- perc_accuracy(train_list, test_list, est_new_features)
 

