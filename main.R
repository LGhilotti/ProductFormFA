######################################################
#### Main script where functions are invoked ##########
#######################################################


####################################################
############## BB with Poisson(lambda) #############
####################################################
set.seed(1234)

# Credible intervals
ci_kmn_poiss_bb <- CI_Kmn_poiss_BB(alpha = - 1, theta = 10, m = 3000, n = 100, 
                                   lambda = 10000, lev = 0.95)
# Plot Credible intervals
plot_Kmn(ci_kmn_poiss_bb)

# Generate from buffet procedure from beginning
# (returns the features, the number of new features for new customers, 
# counts of observed features)
buff_poiss_bb <- buffet_poiss_BB(alpha = - 1, theta = 10, n = 1000, lambda = 1000)

# Matrix of order-of-appearance features from the buffet
ooa_mat_poiss_bb <- create_features_matrix(buff_poiss_bb)

# Generate from buffet procedure given initial n-dimensional sample
# (returns the features and the number of new features for new customers)
buff_poiss_bb_initial_sample <- buffet_poiss_BB_initial_sample(alpha = - 1, theta = 10,
                                                               m=1000, 
                                                               n = length(buff_poiss_bb$num_new),
                                                               counts = buff_poiss_bb$counts, lambda = 1000)

# Empirical Bayes estimate of (alpha, theta, lambda) via EFPF maximization
eb_est_poiss_BB <- EB_EFPF_poiss_BB(n = length(buff_poiss_bb$num_new), 
                        counts = buff_poiss_bb$counts, pars_0 = c(-100, 150, 1000))

eb_est_poiss_BB

###############################################################
################ BB with Negative-Binomial(n*,p) #############
##############################################################
set.seed(1234)

# Credible intervals
ci_kmn_negbin_bb <- CI_Kmn_negbin_BB(alpha = -1, theta = 10, m = 3000, n = 100,
                                     Kn = 10, nstar = 100, p = 0.5, lev = 0.95)
# Plot Credible intervals
plot_Kmn(ci_kmn_negbin_bb)

# Generate from buffet procedure from beginning
# (returns the features and the number of new features for new customers)
buff_negbin_bb <- buffet_negbin_BB(alpha = -5, theta = 10, n = 10000, nstar = 10000, p = 0.5)

# Matrix of order-of-appearance features from the buffet
ooa_mat_negbin_bb <- create_features_matrix(buff_negbin_bb)

# Generate from buffet procedure given initial n-dimensional sample
# (returns the features and the number of new features for new customers)
buff_negbin_bb_initial_sample <- buffet_negbin_BB_initial_sample(alpha=-1,theta=10,m=50, 
                                                                 n= length(buff_negbin_bb$num_new), 
                                                                 counts = buff_negbin_bb$counts, nstar=100, p=0.5)


# Empirical Bayes estimate of (alpha, theta, n*, p) via EFPF maximization
eb_est_negbin_BB <- EB_EFPF_negbin_BB(n = length(buff_negbin_bb$num_new), 
                        counts = buff_negbin_bb$counts, pars_0 = c(-100, 150, 100, 0.5))

eb_est_negbin_BB

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
buff_gamma_ibp <- buffet_gamma_IBP(alpha = 0.6, theta = 10, n = 1000, a = 10, b = 1)

# Matrix of order-of-appearance features from the buffet
ooa_mat_gamma_ibp <- create_features_matrix(buff_gamma_ibp)

# Generate from buffet procedure given initial n-dimensional sample
# (returns the features and the number of new features for new customers)
buff_gamma_ibp_initial_sample <- buffet_gamma_IBP_initial_sample(alpha = - 1, theta = 10,
                                                                 m=100, 
                                                                 n = length(buff_gamma_ibp$num_new),
                                                                 counts = buff_gamma_ibp$counts,
                                                                 a = 1, b= 1)

# Empirical Bayes estimate of (alpha, theta, a, b) via EFPF maximization
eb_est_gamma_IBP <- EB_EFPF_gamma_IBP(n = length(buff_gamma_ibp$num_new), 
                                      counts = buff_gamma_ibp$counts, pars_0 = c(0.5, 50, 10, 10))

eb_est_gamma_IBP
