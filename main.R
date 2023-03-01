######################################################
#### Main script where functions are invoked ##########
#######################################################


## BB with Poisson(lambda)
set.seed(1234)

ci_kmn_poiss_bb <- CI_Kmn_poiss_BB(alpha = - 1, theta = 10, m = 3000, n = 100, 
                                   lambda = 10000, lev = 0.95)
plot_Kmn(ci_kmn_poiss_bb)

buff_poiss_bb <- buffet_poiss_BB(alpha = - 1, theta = 10, n = 100, lambda = 10000)

## BB with Negative-Binomial(n*,p)
set.seed(1234)


ci_kmn_negbin_bb <- CI_Kmn_negbin_BB(alpha = -1, theta = 10, m = 3000, n = 100,
                                     Kn = 10, nstar = 100, p = 0.5, lev = 0.95)
plot_Kmn(ci_kmn_negbin_bb)

buff_negbin_bb <- buffet_negbin_BB(alpha = -1, theta = 10, n = 100, nstar = 100, p = 0.5)

## IBP with Gamma(a,b)
set.seed(1234)
alpha <- 0.2
theta <- 5
a <- 10
b <- 3
lev <- 0.95
n <- 10
Kn <- 3
m <- 10

ci_kmn_gamma_ibp <- CI_Kmn_gamma_IBP(alpha = 0.2, theta = 10, m = 3000, n = 100,
                                     Kn = 10, a = 1, b = 1, lev = 0.95)
plot_Kmn(ci_kmn_gamma_ibp)

buff_gamma_ibp <- buffet_gamma_IBP(alpha = 0.2, theta = 10, n = 100, a = 1, b = 1)



