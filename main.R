######################################################
#### Main script where functions are invoked ##########
#######################################################


####################################################
############## BB with Poisson(lambda) #############
####################################################
set.seed(1234)

# Credible intervals
ci_kmn_poiss_bb <- CI_Kmn_poiss_BB(alpha = - 1, theta = 10, m = 3000, n = 100, 
                                   lambda = 1000, lev = 0.95)
# Plot Credible intervals
plot_Kmn(ci_kmn_poiss_bb)

# Generate from buffet procedure from beginning
# (returns the features, the number of new features for new customers, 
# counts of observed features)
buff_poiss_bb <- buffet_poiss_BB(alpha = - 100, theta = 101, n = 100000, lambda = 100)

# Matrix of order-of-appearance features from the buffet
ooa_mat_poiss_bb <- create_features_matrix(buff_poiss_bb$features)
plot_binary_matrix(ooa_mat_poiss_bb)
plot_binary_matrix(ooa_mat_poiss_bb, max_f = 20)

# Generate from buffet procedure given initial n-dimensional sample
# (returns the features and the number of new features for new customers)
buff_poiss_bb_initial_sample <- buffet_poiss_BB_initial_sample(alpha = - 1, theta = 10,
                                                               m=1000, 
                                                               n = length(buff_poiss_bb$num_new),
                                                               counts = buff_poiss_bb$counts, lambda = 1000)

# Empirical Bayes estimate of (alpha, theta, lambda) via EFPF maximization
eb_efpf_poiss_BB <- EB_EFPF_poiss_BB(n = length(buff_poiss_bb$num_new), 
                        counts = buff_poiss_bb$counts, pars_0 = c(-1, 10, 10))

eb_efpf_poiss_BB

# Empirical Bayes estimate of (alpha, theta, lambda) via MM 
#eb_mm_poiss_BB <- EB_MM_poiss_BB(n = length(buff_poiss_bb$num_new), 
#                                 ntrain = round(length(buff_poiss_bb$num_new)*2/3),
#                                 num_new = buff_poiss_bb$num_new, pars_0 = eb_efpf_poiss_BB)
#
#eb_mm_poiss_BB

n = length(buff_poiss_bb$num_new)
ntrain = round(length(buff_poiss_bb$num_new)*2/3)
ntest = n - ntrain
num_new_m = buff_poiss_bb$num_new[(ntrain+1):n]

alpha = eb_mm_poiss_BB[1]; theta = eb_mm_poiss_BB[2]; lambda = eb_mm_poiss_BB[3];
s=theta+alpha
mm_obj_poiss_BB_rep(pars=c(alpha,s,lambda), ntest=ntest,
                    ntrain = ntrain,
                    num_new_m = num_new_m)

alpha = eb_efpf_poiss_BB[1]; theta = eb_efpf_poiss_BB[2]; lambda = eb_efpf_poiss_BB[3];
s=theta+alpha
mm_obj_poiss_BB_rep(pars=c(alpha,s,lambda), ntest=ntest,
                    ntrain = ntrain,
                    num_new_m = num_new_m)

alpha = -100; theta = 150; lambda = 100; s = alpha + theta
mm_obj_poiss_BB_rep(pars=c(alpha,s,lambda), ntest=ntest,
                    ntrain = ntrain,
                    num_new_m = num_new_m)

alpha = -100; theta = 150; lambda = 10000; s = alpha + theta
mm_obj_poiss_BB_rep(pars=c(alpha,s,lambda), ntest=ntest,
                    ntrain = ntrain,
                    num_new_m = num_new_m)

alpha = -10; theta = 100; lambda = 100; s = alpha + theta
mm_obj_poiss_BB_rep(pars=c(alpha,s,lambda), ntest=ntest,
                    ntrain = ntrain,
                    num_new_m = num_new_m)

alpha = -10; theta = 100; lambda = 10000; s = alpha + theta
mm_obj_poiss_BB_rep(pars=c(alpha,s,lambda), ntest=ntest,
                    ntrain = ntrain,
                    num_new_m = num_new_m)

alpha = -1; theta = 10; lambda = 100; s = alpha + theta
mm_obj_poiss_BB_rep(pars=c(alpha,s,lambda), ntest=ntest,
                    ntrain = ntrain,
                    num_new_m = num_new_m)

alpha = -1; theta = 10; lambda = 10000; s = alpha + theta
mm_obj_poiss_BB_rep(pars=c(alpha,s,lambda), ntest=ntest,
                    ntrain = ntrain,
                    num_new_m = num_new_m)

###############################################################
################ BB with Negative-Binomial(n*,p) #############
##############################################################
set.seed(1234)

# Credible intervals
ci_kmn_negbin_bb <- CI_Kmn_negbin_BB(alpha = -1, theta = 10, m = 3000, n = 100,
                                     Kn = 50, nstar = 1000, p = 0.5, lev = 0.95)
# Plot Credible intervals
plot_Kmn(ci_kmn_negbin_bb)

# Generate from buffet procedure from beginning
# (returns the features and the number of new features for new customers)
buff_negbin_bb <- buffet_negbin_BB(alpha = -1, theta = 10, n = 100000, nstar = 1000, p = 0.5)
n_feat <- length(buff_negbin_bb$counts)

# Matrix of order-of-appearance features from the buffet
ooa_mat_negbin_bb <- create_features_matrix(buff_negbin_bb$features)
plot_binary_matrix(ooa_mat_negbin_bb)
plot_binary_matrix(ooa_mat_negbin_bb, max_f = 20)


# Generate from buffet procedure given initial n-dimensional sample
# (returns the features and the number of new features for new customers)
buff_negbin_bb_initial_sample <- buffet_negbin_BB_initial_sample(alpha=-1,theta=10,m=50, 
                                                                 n= length(buff_negbin_bb$num_new), 
                                                                 counts = buff_negbin_bb$counts, nstar=100, p=0.5)


# Empirical Bayes estimate of (alpha, theta, n*, p) via EFPF maximization
eb_efpf_negbin_BB <- EB_EFPF_fixed_negbin_BB(n = length(buff_negbin_bb$num_new), 
                                        counts = buff_negbin_bb$counts, 
                                        pars_0 = c(-10, 15, 50, 0.2),
                                        opt = c(T,T,T,T))

# Empirical Bayes estimate of (alpha, theta, mu, sigma2) via EFPF maximization (m-v param)
eb_efpf_negbin_mv_BB <- EB_EFPF_fixed_negbin_mv_BB(n = length(buff_negbin_bb$num_new), 
                                              counts = buff_negbin_bb$counts, 
                                              pars_0 = c(-1, 15, 10, true_var),
                                              opt = c(T,T,T,F))

########################################################################################
############# START CHECK ##############################################################
########################################################################################

# Empirical Bayes estimate of (alpha, theta, mu, sigma2) via EFPF maximization (m-v param)
# with replicates (JUST FOR A CHECK)
set.seed(1234)
n_replicates <- 100
l_counts_negbin_bb <- vector("list",n_replicates) # list of counts of features for each replicate
alpha = -100; theta = 101; n = 10000; nstar = 100; p = 0.8

true_mean <- nstar*(1-p)/p
true_mean

true_var <- nstar*(1-p)/(p**2)
true_var

# Generate from buffet procedure multiple times and store in list
for (i in 1:n_replicates){
  
  buff_negbin_bb_rep <- buffet_negbin_BB(alpha, theta, n, nstar, p)
  l_counts_negbin_bb[[i]] <- buff_negbin_bb_rep$counts
}

avg_n_obs_feat <- length(unlist(l_counts_negbin_bb))/n_replicates
avg_n_obs_feat
eb_efpf_mv_replicates <- EB_EFPF_fixed_negbin_mv_BB_replicates(n, 
                                                               counts = l_counts_negbin_bb,
                                                               pars_0 = c(-1, 5, 10, 20),
                                                               opt = c(T,T,T,T))
eb_efpf_mv_replicates

# -> single realization of the sample
eb_efpf_mv_single <- EB_EFPF_fixed_negbin_mv_BB(n,
                                                counts = l_counts_negbin_bb[[1]],
                                                pars_0 = c(-3, 10, 30, 60),
                                                opt = c(T,T,T,T))
eb_efpf_mv_single
#################################################################################
############### END CHECK #######################################################
#################################################################################

###############################################################
############ IBP with Gamma(a,b) ##############################
###############################################################
set.seed(1234)

# Credible intervals
ci_kmn_gamma_ibp <- CI_Kmn_gamma_IBP(alpha = 0.6, theta = 10, m = 3000, n = 100,
                                     Kn = 10, a = 10, b = 1, lev = 0.95)
# Plot Credible intervals
plot_Kmn(ci_kmn_gamma_ibp)

# Generate from buffet procedure from beginning
# (returns the features and the number of new features for new customers)
buff_gamma_ibp <- buffet_gamma_IBP(alpha = 0.2, theta = 2, n = 10000, a = 2, b = 1)

# Matrix of order-of-appearance features from the buffet
ooa_mat_gamma_ibp <- create_features_matrix(buff_gamma_ibp$features)
plot_binary_matrix(ooa_mat_gamma_ibp)
plot_binary_matrix(ooa_mat_gamma_ibp, max_f = 20)

# Generate from buffet procedure given initial n-dimensional sample
# (returns the features and the number of new features for new customers)
buff_gamma_ibp_initial_sample <- buffet_gamma_IBP_initial_sample(alpha = - 1, theta = 10,
                                                                 m=100, 
                                                                 n = length(buff_gamma_ibp$num_new),
                                                                 counts = buff_gamma_ibp$counts,
                                                                 a = 1, b= 1)

# Empirical Bayes estimate of (alpha, theta, a, b) via EFPF maximization
eb_efpf_gamma_IBP <- EB_EFPF_gamma_IBP(n = length(buff_gamma_ibp$num_new), 
                                      counts = buff_gamma_ibp$counts, pars_0 = c(0.5, 10, 5, 5))

eb_efpf_gamma_IBP



