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
# (returns the features and the number of new features for new customers)
buff_poiss_bb <- buffet_poiss_BB(alpha = - 1, theta = 10, n = 10, lambda = 1000)

# Matrix of order-of-appearance features from the buffet
ooa_mat_poiss_bb <- create_features_matrix(buff_poiss_bb)

# Generate from buffet procedure given initial n-dimensional sample
# (returns the features and the number of new features for new customers)
counts <- c(10,8,6,5,4,5,3,2,3,3)
buff_poiss_bb_initial_sample <- buffet_poiss_BB_initial_sample(alpha = - 1, theta = 10,
                                                               m=100, n = 10,
                                                               counts = counts, lambda = 1000)

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
buff_negbin_bb <- buffet_negbin_BB(alpha = -1, theta = 10, n = 100, nstar = 100, p = 0.5)

# Matrix of order-of-appearance features from the buffet
ooa_mat_negbin_bb <- create_features_matrix(buff_negbin_bb)

# Generate from buffet procedure given initial n-dimensional sample
# (returns the features and the number of new features for new customers)
counts <- c(10,8,6,5,4,5,3,2,3,3)
buff_negbin_bb_initial_sample <- buffet_negbin_BB_initial_sample(alpha = - 1, theta = 10,
                                                               m=100, n = 10,
                                                               counts = counts,
                                                               nstar = 100, p=0.5)



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
buff_gamma_ibp <- buffet_gamma_IBP(alpha = 0.6, theta = 10, n = 100, a = 10, b = 1)

# Matrix of order-of-appearance features from the buffet
ooa_mat_gamma_ibp <- create_features_matrix(buff_gamma_ibp)

# Generate from buffet procedure given initial n-dimensional sample
# (returns the features and the number of new features for new customers)
counts <- c(10,8,6,5,4,5,3,2,3,3)
buff_gamma_ibp_initial_sample <- buffet_gamma_IBP_initial_sample(alpha = - 1, theta = 10,
                                                                 m=100, n = 10,
                                                                 counts = counts,
                                                                 a = 1, b= 1)


