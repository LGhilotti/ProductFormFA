#
#### How to use the functions? ------
#
rm(list = ls())
library(ProductFormFA)

source("R_script_paper/utils.R", echo=TRUE)

# 1) Random generation of the processes ----
set.seed(12345)
seed <- 12345

n <- 500
alpha_true <- -1
theta_true <- 10
N_true <- 30

# 1.1) Beta-Bernoulli with Poisson(lambda) mixture
# sample_PoissonBB <- rPoissonBB(alpha = alpha_true, theta = theta_true, lambda = N_true, n = 1000)
# feature_matrix_PoissonBB <- convert_features_matrix(sample_PoissonBB$features)



set.seed(012345)
probs1 <- rbeta(N_true, -alpha_true, alpha_true + theta_true)
probs1 <- c(t(replicate(n = n, expr = probs1)))

hist(probs1)

probs2 <- sample(c(0.1, 0.5, 0.25), N_true, replace = TRUE)
probs2 <- c(t(replicate(n = n, expr = probs2)))

hist(probs2)

feature_matrix_PoissonBB <- cbind(
  matrix(rbinom(n * N_true, 1, prob = probs2), nrow = n, ncol = N_true)  
)

feature_matrix_PoissonBB <- cbind(matrix(rbinom(50 * 10, size = 1, prob = 0.9), 50, 20),
                       matrix(rbinom(50 * 10, size = 1, prob = 0.9), 50, 20),
                       matrix(rbinom(50 * 10, size = 1, prob = 0.9), 50, 20),
                       matrix(rbinom(50 * 10, size = 1, prob = 0.8), 50, 20),
                       matrix(rbinom(50 * 10, size = 1, prob = 0.7), 50, 20),
                       matrix(rbinom(50 * 10, size = 1, prob = 0.7), 50, 20),
                       matrix(rbinom(50 * 10, size = 1, prob = 0.01), 50, 20),
                       matrix(rbinom(50 * 10, size = 1, prob = 0.01), 50, 20),
                       matrix(rbinom(50 * 100, size = 1, prob = 0.5), 50, 100))

feature_matrix_PoissonBB <- BCI > 0#feature_matrix_PoissonBB[, colSums(feature_matrix_PoissonBB) > 0]

plot(rarefaction(feature_matrix_PoissonBB))

hist(colMeans(feature_matrix_PoissonBB))

plot(table(colSums(feature_matrix_PoissonBB)))

# -\alpha ~ Gamma(a_alpha, b_alpha)
# s ~ Gamma(a_s, b_s)
# s = \alpha + \theta

# Hyperparameters
hyper <- list(
  a_alpha = 0.01, b_alpha = 0.01,
  a_s = 0.01, b_s = 0.01,
  lambda = beta_binomial_estimator(feature_matrix_PoissonBB)
)

prior_obj <- prior(model = "PoissonBB", hyper = hyper)
summary(prior_obj)

# Initialization
init <- list(alpha_0 = -1, s_0 = 1)
init_obj <- initialization(model = "PoissonBB", init = init)

# MCMC setting
# tau = "tuning parameter for MCMC"
mcmcparams <- list(tau = 0.05, S = 5000, n_burnin = 1000, thin = 5)
mcmcparams_obj <- mcmcparameters(model = "PoissonBB", mcmcparams = mcmcparams)

# Fit the model
PoissonBB_fit <- GibbsFA(
  feature_matrix = feature_matrix_PoissonBB,
  model = "PoissonBB",
  prior = prior_obj,
  initialization = init_obj,
  mcmcparams = mcmcparams_obj
)

# 2.1.B) Beta-Bernoulli with Poisson(lambda) mixture - EB version

# EB parameters
eb_init <- list(alpha = -1, s = 1)
eb_known <- list(lambda = beta_binomial_estimator(feature_matrix_PoissonBB))
eb_params_obj <- eb_params(model = "PoissonBB", init = eb_init, known = eb_known)

# Fit the model with EB
eb_PoissonBB_fit <- GibbsFA_eb(
  feature_matrix = feature_matrix_PoissonBB,
  model = "PoissonBB",
  eb_params = eb_params_obj
)


data.frame(
  truth = c(alpha_true, theta_true, N_true, -alpha_true / theta_true),
  EB = c(eb_PoissonBB_fit$alpha, eb_PoissonBB_fit$theta, eb_PoissonBB_fit$lambda, -eb_PoissonBB_fit$alpha / eb_PoissonBB_fit$theta),
  Bayes = c(mean(PoissonBB_fit$alpha_chain), mean(PoissonBB_fit$theta_chain), eb_PoissonBB_fit$lambda, mean(-PoissonBB_fit$alpha_chain / PoissonBB_fit$theta_chain))
)

plot(PoissonBB_fit$alpha_chain, type = "l")
hist(PoissonBB_fit$alpha_chain, breaks = 100)
abline(v = eb_PoissonBB_fit$alpha, lty = "dotted")
abline(v = , lty = "dashed", col = "red")

hist(PoissonBB_fit$theta_chain, breaks = 100)
abline(v = eb_PoissonBB_fit$theta, lty = "dotted")
abline(v = mean(PoissonBB_fit$theta_chain), lty = "dashed", col = "red")

hist(-PoissonBB_fit$alpha_chain / PoissonBB_fit$theta_chain, breaks = 100)

plot(x = PoissonBB_fit, type = "rarefaction", n_reorderings = 10, seed = seed) # option: seed
plot(x = eb_PoissonBB_fit, type = "rarefaction", n_reorderings = 10) # option: seed

plot(x = eb_PoissonBB_fit, type = "extrapolation", M = 200, n_reorderings = 10) # option: seed


#
# # 2.2) Beta-Bernoulli with NB(n0, mu0) mixture
#
# feature_matrix_NegBinBB <- convert_features_matrix(sample_NegBinBB$features)
#
# # Hyperparameters
# hyper <- list(
#   a_alpha = 1, b_alpha = 0.1,
#   a_s = 2, b_s = 0.2,
#   n0 = 10, mu0 = 10
# )
#
# prior_obj <- prior(model = "NegBinBB", hyper = hyper)
#
# summary(prior_obj)
#
# # Initialization
# init <- list(alpha_0 = -1, s_0 = 1)
#
# init_obj <- initialization(model = "NegBinBB", init = init)
#
# # MCMC setting
# mcmcparams <- list(tau = 0.1, S = 100, n_burnin = 10, thin = 2)
#
# mcmcparams_obj <- mcmcparameters(model = "NegBinBB", mcmcparams = mcmcparams)
#
# # Fit the model
# NegBinBB_fit <- GibbsFA(
#   feature_matrix = feature_matrix_NegBinBB,
#   model = "NegBinBB",
#   prior = prior_obj,
#   initialization = init_obj,
#   mcmcparams = mcmcparams_obj
# )
#
#
# # 2.2.B1) Beta-Bernoulli with NB(n0,mu0) mixture - EB version
#
# # EB parameters
# eb_init <- list(alpha = -1, s = 1, mu0 = 10)
# eb_known <- list(n0 = 1000)
# eb_params_obj <- eb_params(model = "NegBinBB", init = eb_init, known = eb_known)
#
# # Fit the model with EB
# eb_NegBinBB_fit <- GibbsFA_eb(
#   feature_matrix = feature_matrix_NegBinBB,
#   model = "NegBinBB",
#   eb_params = eb_params_obj
# )
#
#
# # 2.2.B2) Beta-Bernoulli with NB(n0,p) mixture - EB version
#
# # EB parameters
# eb_init <- list(alpha = -1, s = 1, p = 0.8)
# eb_known <- list(n0 = 10)
# eb_params_obj <- eb_params(model = "NegBinBB_np", init = eb_init, known = eb_known)
#
# # Fit the model with EB
# eb_NegBinBB_np_fit <- GibbsFA_eb(
#   feature_matrix = feature_matrix_NegBinBB,
#   model = "NegBinBB_np",
#   eb_params = eb_params_obj
# )
#

# 2.3) IBP with Gamma(a,b) mixture

feature_matrix_GammaIBP <- convert_features_matrix(sample_GammaIBP$features)

# Hyperparameters
hyper <- list(
  a_alpha = 1, b_alpha = 2,
  a_s = 2, b_s = 0.2,
  q = 0.5, r = 1, t = 1
)

prior_obj <- prior(model = "GammaIBP", hyper = hyper)

summary(prior_obj)

# Initialization
init <- list(alpha_0 = 0.5, s_0 = 1, a_0 = 1, b_0 = 1)

init_obj <- initialization(model = "GammaIBP", init = init)

# MCMC setting
mcmcparams <- list(
  sigq_alpha = 0.1, sigq_s = 1,
  S = 10000, n_burnin = 1000, thin = 2
)

mcmcparams_obj <- mcmcparameters(model = "GammaIBP", mcmcparams = mcmcparams)

# Fit the model
GammaIBP_fit <- GibbsFA(
  feature_matrix = feature_matrix_PoissonBB,
  model = "GammaIBP",
  prior = prior_obj,
  initialization = init_obj,
  mcmcparams = mcmcparams_obj
)


# 2.3.B) IBP with Gamma(a,b) mixture - EB version

# EB parameters
eb_init <- list(alpha = 0.5, s = 1, a = 10, b = 0.5)
eb_known <- list()
eb_params_obj <- eb_params(model = "GammaIBP", init = eb_init, known = eb_known)

# Fit the model with EB
eb_GammaIBP_fit <- GibbsFA_eb(
  feature_matrix = feature_matrix_PoissonBB,
  model = "GammaIBP",
  eb_params = eb_params_obj
)

#
#
# # 3) Estimate the total richness ----
#
#
#
# # 3.2) Beta-Bernoulli with NB(n0,mu0) mixture
#
# NegBinBB_rich <- total_richness(object = NegBinBB_fit)
#
# # 3.2.B) Beta-Bernoulli with NB(n0,mu0) mixture - EB version
#
# eb_NegBinBB_rich <- total_richness(object = eb_NegBinBB_fit)
# # these are the NB parameters of N' (in n0,mu0 parametrization), then add Kn
#
#
#
# # 4) Rarefaction ----
#
# # 4.1) Beta-Bernoulli with Poisson(lambda) mixture
#
# PoissonBB_rare <- rarefaction(object = PoissonBB_fit) # option: seed
#
# # 4.1.B) Beta-Bernoulli with Poisson(lambda) mixture - EB version
#
# eb_PoissonBB_rare <- rarefaction(object = eb_PoissonBB_fit)
#
# # 4.2) Beta-Bernoulli with NB(n0,mu0) mixture
#
# NegBinBB_rare <- rarefaction(object = NegBinBB_fit) # option: seed
#
# # 4.2.B) Beta-Bernoulli with NB(n0,mu0) mixture - EB version
#
# eb_NegBinBB_rare <- rarefaction(object = eb_NegBinBB_fit) #
#
# # 4.3) IBP with Gamma(a,b) mixture
#
# GammaIBP_rare <- rarefaction(object = GammaIBP_fit) # option: seed
#
# # 4.3.B) IBP with Gamma(a,b) mixture - EB version
#
# eb_GammaIBP_rare <- rarefaction(object = eb_GammaIBP_fit) # option: seed
#
#
#
# # 5) Extrapolation ----
#
# # 5.1) Beta-Bernoulli with Poisson(lambda) mixture
#
# PoissonBB_extr <- extrapolation(object = PoissonBB_fit, M = 20) # option: only_last = TRUE, seed
#
# # 5.1.B) Beta-Bernoulli with Poisson(lambda) mixture - EB version
#
# eb_PoissonBB_extr <- extrapolation(object = eb_PoissonBB_fit, M = 20) # option: only_last = TRUE
# # this is the Poisson parameter of N', then add Kn
#
# # 5.2) Beta-Bernoulli with NB(n0,mu0) mixture
#
# NegBinBB_extr <- extrapolation(object = NegBinBB_fit, M = 20) # option: only_last = TRUE, seed
#
# # 5.2.B) Beta-Bernoulli with NB(n0,mu0) mixture - EB version
#
# eb_NegBinBB_extr <- extrapolation(object = eb_NegBinBB_fit, M = 20) # option: only_last = TRUE
# # these are the NB parameters of N' (in n0,mu0 parametrization), then add Kn
#
# # 5.3) IBP with Gamma(a,b) mixture
#
# GammaIBP_extr <- extrapolation(object = GammaIBP_fit, M = 20) # option: only_last = TRUE, seed
#
# # 5.3.B) IBP with Gamma(a,b) mixture - EB version
#
# eb_GammaIBP_extr <- extrapolation(object = eb_GammaIBP_fit, M = 20) # option: only_last = TRUE
# # these are the NB parameters of N' (in n0,mu0 parametrization), then add Kn
#
#
#
#
# # 6) Plot richness ----
#
# # 6.1) Beta-Bernoulli with Poisson(lambda) mixture
#
# plot(x = PoissonBB_fit, type = "richness") # option: bw (bandwidth for density estimation)
#
# # 6.1.B) Beta-Bernoulli with Poisson(lambda) mixture - EB version
#
# plot(x = eb_PoissonBB_fit, type = "richness")
#
# # 6.2) Beta-Bernoulli with NB(n0,mu0) mixture
#
# plot(x = NegBinBB_fit, type = "richness") # option: bw (bandwidth for density estimation)
#
# # 6.2.B) Beta-Bernoulli with NB(n0,mu0) mixture - EB version
#
# plot(x = eb_NegBinBB_fit, type = "richness")
#
#
# # 7) Plot Rarefaction ----
#
#
# # 7.2) Beta-Bernoulli with NB(n0,mu0) mixture
#
# plot(x = NegBinBB_fit, type = "rarefaction", n_reorderings = 10) # option: seed
#
# # 7.2.B) Beta-Bernoulli with NB(n0,mu0) mixture - EB version
#
# plot(x = eb_NegBinBB_fit, type = "rarefaction", n_reorderings = 10) # option: seed
#
# 7.3) IBP with Gamma(a,b) mixture

plot(x = GammaIBP_fit, type = "rarefaction", n_reorderings = 10) # option: seed

# 7.3.B) IBP with Gamma(a,b) mixture - EB version

plot(x = eb_GammaIBP_fit, type = "rarefaction", n_reorderings = 10) # option: seed

#
#
# # 8) Plot Extrapolation ----
#
# # 8.1) Beta-Bernoulli with Poisson(lambda) mixture
#
# plot(x = PoissonBB_fit, type = "extrapolation", M = 20, n_reorderings = 10) # option: seed
#
# # 8.1.B) Beta-Bernoulli with Poisson(lambda) mixture - EB version
#
# plot(x = eb_PoissonBB_fit, type = "extrapolation", M = 20, n_reorderings = 10) # option: seed
#
# # 8.2) Beta-Bernoulli with NB(n0,mu0) mixture
#
# plot(x = NegBinBB_fit, type = "extrapolation", M = 20, n_reorderings = 10) # option:  seed
#
# # 8.2.B) Beta-Bernoulli with NB(n0,mu0) mixture - EB version
#
# plot(x = eb_NegBinBB_fit, type = "extrapolation", M = 20, n_reorderings = 10) # option: seed
#
# 8.3) IBP with Gamma(a,b) mixture

plot(x = GammaIBP_fit, type = "extrapolation", M = 200, n_reorderings = 10) # option:  seed

# 8.3.B) IBP with Gamma(a,b) mixture - EB version

plot(x = eb_GammaIBP_fit, type = "extrapolation", M = 20, n_reorderings = 10) # option: seed

data(dune.env)
pool <- with(dune.env, specpool(dune, Management))
pool
