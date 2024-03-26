#
#### How to use the functions? ------
#

library(ProductFormFA)

# 1) Random generation of the processes ----

seed = 1234

# 1.1) Beta-Bernoulli with Poisson(lambda) mixture

sample_PoissonBB <- rPoissonBB(alpha = -1, theta = 2, lambda = 10, n = 10)

sample_PoissonBB

# 1.2) Beta-Bernoulli with NB(nstar, p) mixture

sample_NegBinBB <- rNegBinBB(alpha = -1, theta = 2, nstar = 10, p = 0.5, n = 10)

sample_NegBinBB

# 1.3) IBP with Gamma(a, b) mixture

sample_GammaIBP <- rGammaIBP(alpha = 0.5, theta = 2, a = 1, b = 1, n = 10)

sample_GammaIBP


