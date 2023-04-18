##############################################################
############ REAL DATA APPLICATION: FUNGI ###################
#############################################################

data_fungi_all <- read.csv("data/data_fungi/d_fungi.csv")
data_fungi <- subset(data_fungi_all, select = - c(LogID,DC,readcount,unk))

# number of observed fungi species
n_observed_species <- sum(colSums(data_fungi) !=0)

# binary matrix
fungi_mat <- as.matrix((data_fungi!=0)+0)
rowSums(fungi_mat)

# number of sites
L = nrow(fungi_mat)

fungi_list <- create_features_list(fungi_mat)
plot_trajectory(fungi_list)

# set the dimension of the training set and the training set itself
Ns <- c(10, 20, 50)
# list of ggplots - poiss
gg_list_poiss <- vector("list", length = length(Ns) )
# list of ggplots - negbin
gg_list_negbin <- vector("list", length = length(Ns) )

for (j in 1:length(Ns)){
  N <- Ns[j]
  
  train_list <- fungi_list[1:N]
  train_mat <- fungi_mat[1:N, ]
  
  # set the dimension of the test set and the test set itself
  M <- L - N
  test_list <- fungi_list[(N+1):L]
  test_mat <- fungi_mat[(N+1):L, ]
  
  # compute number of observed feature in the training set Kn
  Kn <- length(unique(unlist(train_list)))
  
  # EB: fix the model hyperparameters based on EFPF-maximization
  train_counts <- colSums(train_mat)[colSums(train_mat)!=0]
  
  # Fit Poisson model
  # EB: fix the model hyperparameters based on EFPF-maximization
  eb_efpf_poiss_BB <- EB_EFPF_poiss_BB(n = N, 
                                       counts = train_counts, pars_0 = c(-2, 5, 10))
  alpha_est_poiss_efpf <- eb_efpf_poiss_BB[1]
  theta_est_poiss_efpf <- eb_efpf_poiss_BB[2]
  lambda_est_poiss_efpf <- eb_efpf_poiss_BB[3]
  
  # Fit Negative-Binomial model
  # EB: fix the model hyperparameters based on EFPF-maximization
  eb_efpf_negbin_BB <- EB_EFPF_fixed_negbin_mv_BB(n = N, 
                                                  counts = train_counts, pars_0 = c(-2, 5, 100, 150))
  
  alpha_est_negbin_efpf <- eb_efpf_negbin_BB[1]
  theta_est_negbin_efpf <- eb_efpf_negbin_BB[2]
  mu_est_negbin_efpf <- eb_efpf_negbin_BB[3]
  sigma2_est_negbin_efpf <- eb_efpf_negbin_BB[4]
  
  nstar_est_negbin_efpf <- (mu_est_negbin_efpf**2)/(sigma2_est_negbin_efpf - mu_est_negbin_efpf)
  p_est_negbin_efpf <- mu_est_negbin_efpf / sigma2_est_negbin_efpf
  
  # Compute E[KMN|Z_N] together with 90% credible intervals, for given N and train set
  CI_given_sample_poiss_efpf <- CI_Kmn_poiss_BB(alpha_est_poiss_efpf, theta_est_poiss_efpf,
                                                M, N, lambda_est_poiss_efpf, 0.95)
  
  CI_given_sample_negbin_efpf <- CI_Kmn_negbin_BB(alpha_est_negbin_efpf, theta_est_negbin_efpf,
                                                  M, N, Kn, nstar_est_negbin_efpf, 
                                                  p_est_negbin_efpf, 0.95)
  
  gg_list_poiss[[j]] <- plot_Kmn_given_sample_with_observed(N, fungi_list, CI_given_sample_poiss_efpf)

  gg_list_negbin[[j]] <- plot_Kmn_given_sample_with_observed(N, fungi_list, CI_given_sample_negbin_efpf)
  
}

# display plots on the same page
ggarrange(gg_list_poiss[[1]], gg_list_poiss[[2]], gg_list_poiss[[3]],
          gg_list_negbin[[1]], gg_list_negbin[[2]], gg_list_negbin[[3]],
          ncol = 3, nrow = 2)

