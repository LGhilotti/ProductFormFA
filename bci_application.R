##############################################################
############ REAL DATA APPLICATION: BCI (vegan) ###################
#############################################################

library(vegan)
data(BCI)
bci_mat <- as.matrix((BCI > 0) +0)
freq <- colSums(bci_mat)
# to check that all species in the dataset are actually present 
sum(freq ==0) # must be =0

# number of sites
L = nrow(bci_mat)

bci_list <- create_features_list(bci_mat)
plot_trajectory(bci_list)

# set the dimension of the training set (just to plot)
Ns <- c(10, 20, 30)
# list of ggplots - poiss
gg_list_poiss <- vector("list", length = length(Ns) )
# list of ggplots - negbin
gg_list_negbin <- vector("list", length = length(Ns) )

# EB: fix the model hyperparameters based on EFPF-maximization
bci_counts <- colSums(bci_mat)[colSums(bci_mat)!=0]

set.seed(1234)

# Fit Poisson model
# EB: fix the model hyperparameters based on EFPF-maximization
eb_efpf_poiss_BB <- EB_EFPF_poiss_BB(n = L, 
                                     counts = bci_counts, pars_0 = c(-2, 5, 10))
alpha_est_poiss_efpf <- eb_efpf_poiss_BB[1]
theta_est_poiss_efpf <- eb_efpf_poiss_BB[2]
lambda_est_poiss_efpf <- eb_efpf_poiss_BB[3]

# Fit Negative-Binomial model
# EB: fix the model hyperparameters based on EFPF-maximization
eb_efpf_negbin_BB <- EB_EFPF_fixed_negbin_mv_BB(n = L, 
                                                counts = bci_counts, pars_0 = c(-2, 5, 100, 150))

alpha_est_negbin_efpf <- eb_efpf_negbin_BB[1]
theta_est_negbin_efpf <- eb_efpf_negbin_BB[2]
mu_est_negbin_efpf <- eb_efpf_negbin_BB[3]
sigma2_est_negbin_efpf <- eb_efpf_negbin_BB[4]

nstar_est_negbin_efpf <- (mu_est_negbin_efpf**2)/(sigma2_est_negbin_efpf - mu_est_negbin_efpf)
p_est_negbin_efpf <- mu_est_negbin_efpf / sigma2_est_negbin_efpf

for (j in 1:length(Ns)){
  N <- Ns[j]
  
  # set the dimension of the test set and the test set itself
  M <- L - N
  
  # compute number of observed feature in the training set Kn
  Kn <- length(unique(unlist(bci_list[1:N])))
  
  # Compute E[KMN|Z_N] together with 90% credible intervals, for given N and train set
  CI_given_sample_poiss_efpf <- CI_Kmn_poiss_BB(alpha_est_poiss_efpf, theta_est_poiss_efpf,
                                                M, N, lambda_est_poiss_efpf, 0.95)
  
  CI_given_sample_negbin_efpf <- CI_Kmn_negbin_BB(alpha_est_negbin_efpf, theta_est_negbin_efpf,
                                                  M, N, Kn, nstar_est_negbin_efpf, 
                                                  p_est_negbin_efpf, 0.95)
  
  gg_list_poiss[[j]] <- plot_Kmn_given_sample_with_observed(N, bci_list, CI_given_sample_poiss_efpf)
  
  gg_list_negbin[[j]] <- plot_Kmn_given_sample_with_observed(N, bci_list, CI_given_sample_negbin_efpf)
  
}

# display plots on the same page
ggarrange(gg_list_poiss[[1]], gg_list_poiss[[2]], gg_list_poiss[[3]],
          gg_list_negbin[[1]], gg_list_negbin[[2]], gg_list_negbin[[3]],
          ncol = 3, nrow = 2)
