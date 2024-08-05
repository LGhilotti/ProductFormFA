#
#### Example on "how to use the available routines" ------
#

rm(list = ls())
library(ProductFormFA)
library(tidyverse)

# 1) Random generation of the processes ----

seed = 1234

# 1.1) Beta-Bernoulli with Poisson(lambda) mixture

sample_PoissonBB <- rPoissonBB(alpha = -1, theta = 2, lambda = 10, n = 10)

sample_PoissonBB

feature_matrix_PoissonBB <- convert_features_matrix(sample_PoissonBB$features)

# 1.2) Beta-Bernoulli with NB(n0, mu0) mixture

sample_NegBinBB <- rNegBinBB(alpha = -1, theta = 12, n0 = 1000, mu0 = 100, n = 100)

sample_NegBinBB

feature_matrix_NegBinBB <- convert_features_matrix(sample_NegBinBB$features)

# 1.3) IBP with Gamma(a, b) mixture

sample_GammaIBP <- rGammaIBP(alpha = 0.5, theta = 10, a = 1, b = 1, n = 20)

sample_GammaIBP

feature_matrix_GammaIBP <- convert_features_matrix(sample_GammaIBP$features)


# 2) Run the models on generated data ----

## 2A) Empirical Bayes approach ----

### 2.1) Beta-Bernoulli with Poisson(lambda) mixture - EB version

# EB parameters 
eb_init <- list(alpha = -1, s = 1)
eb_known <- list(lambda = 10)
eb_params_obj <- eb_params(model = "PoissonBB", init = eb_init, known = eb_known )

# Fit the model with EB
eb_PoissonBB_fit <- GibbsFA_eb(feature_matrix = feature_matrix_PoissonBB, 
                               model = "PoissonBB", type = "EFPF",
                               eb_params =  eb_params_obj)


### 2.2) Beta-Bernoulli with NB(n0,mu0) mixture - EB version

# EB parameters 
eb_init <- list(alpha = -1, s = 1 )
eb_known <- list( mu0 = 100, var_fct = 100)
eb_params_obj <- eb_params(model = "NegBinBB", init = eb_init, known = eb_known )

# Fit the model with EB
eb_NegBinBB_fit <- GibbsFA_eb(feature_matrix = feature_matrix_NegBinBB, 
                              model = "NegBinBB", type = "EFPF",
                              eb_params =  eb_params_obj )


### 2.3) IBP with Gamma(a,b) mixture - EB version

# EB parameters 
eb_init <- list(alpha = 0.5, s = 1, Gamma = 1 )
eb_known <- list()
eb_params_obj <- eb_params(model = "IBP", init = eb_init, known = eb_known )

# Fit the model with EB
eb_GammaIBP_fit <- GibbsFA_eb(feature_matrix = feature_matrix_GammaIBP, 
                              model = "GammaIBP", type = "EFPF",
                              eb_params =  eb_params_obj,
                              var_GammaIBP = 10)



## 2B) Fully-Bayesian approach ----

### 2.2) Beta-Bernoulli with NB(n0, mu0) mixture

# Hyperparameters
hyper <- list(a_alpha = 0.1, b_alpha = 0.1,
              a_s = 0.1, b_s = 0.1,
              n0 = 10, mu0 = 10)
prior_obj <- prior(model = "NegBinBB", hyper = hyper) 

# Initialization 
init <- list(alpha_0 = -1, s_0 = 1)
init_obj <- initialization(model = "NegBinBB", init = init )

# MCMC setting
mcmcparams <- list(tau = 0.1, S = 1000, n_burnin = 100, thin = 2)
mcmcparams_obj <- mcmcparameters(model = "NegBinBB", mcmcparams = mcmcparams)

# Fit the model
NegBinBB_fit <- GibbsFA(feature_matrix = feature_matrix_NegBinBB, 
                        model = "NegBinBB", 
                        prior = prior_obj, 
                        initialization = init_obj, 
                        mcmcparams = mcmcparams_obj)



### 2.3) IBP with Gamma(a,b) mixture

# Hyperparameters
hyper <- list(a = 1, b = 1,
              a_alpha = 0.1, b_alpha = 0.1,
              a_s = 0.1 , b_s = 0.1)
prior_obj <- prior(model = "GammaIBP_single_prior", hyper = hyper) 

# Initialization 
init <- list(alpha_0 = 0.5, s_0 = 1)
init_obj <- initialization(model = "GammaIBP_single_prior", init = init )

# MCMC setting
mcmcparams <- list(sigq_alpha = 0.1, sigq_s = 1, 
                   S = 1000, n_burnin = 100, thin = 2)
mcmcparams_obj <- mcmcparameters(model = "GammaIBP", mcmcparams = mcmcparams)

# Fit the model
GammaIBP_fit <- GibbsFA(feature_matrix = feature_matrix_GammaIBP, 
                        model = "GammaIBP_single_prior", 
                        prior = prior_obj, 
                        initialization = init_obj, 
                        mcmcparams = mcmcparams_obj)




# 3) Estimate the total richness ----

## 3A) Empirical Bayes approach ----

# 3.1) Beta-Bernoulli with Poisson(lambda) mixture - EB version
Kn_PoissonBB <- ncol(feature_matrix_PoissonBB[, colSums(feature_matrix_PoissonBB)!=0])

params_richness_EFPF_PoissonBB <- tibble( lambda_prime = total_richness(eb_PoissonBB_fit)$lambda_post) %>%
  add_column(Model = "Poisson BB") %>%
  mutate(lb = Kn_PoissonBB + qpois(0.025, lambda_prime, lower.tail = TRUE, log.p = FALSE),
         ub = Kn_PoissonBB + qpois(0.975, lambda_prime, lower.tail = TRUE, log.p = FALSE),
         mean = Kn_PoissonBB + lambda_prime) %>%
  select(mean, lb, ub, Model)


# 3.2) Beta-Bernoulli with NB(n0,mu0) mixture - EB version
Kn_NegBinBB <- ncol(feature_matrix_NegBinBB[, colSums(feature_matrix_NegBinBB)!=0])

params_richness_EFPF_NegBinBB <- tibble( n0_prime = 
                                               total_richness(eb_NegBinBB_fit)$n0_post,
                                             mu0_prime = 
                                               total_richness(eb_NegBinBB_fit)$mu0_post) %>%
  add_column(Model = paste0("NegBinomial BB x", eb_NegBinBB_fit$var_fct)) %>%
  mutate(p_prime = 1/(mu0_prime/n0_prime + 1),
         lb = Kn_NegBinBB + qnbinom(0.025, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE),
         ub = Kn_NegBinBB + qnbinom(0.975, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE),
         mean = Kn_NegBinBB + mu0_prime ) %>%
  select(mean, lb, ub, Model)


## 3B) Fully-Bayesian approach ----

# 3.2) Beta-Bernoulli with NB(n0,mu0) mixture
samples_richness <- total_richness(NegBinBB_fit)

params_richness_EFPF_NegBinBB <- tibble( mean = mean(samples_richness),
                                         lb = quantile(samples_richness, probs = 0.025 ),
                                         ub = quantile(samples_richness, probs = 0.975 )) %>%
  add_column(Model = paste0("NegBinomial BB x", eb_NegBinBB_fit$var_fct)) %>%
  select(mean, lb, ub, Model)



# 4) Rarefaction ----

# 4.1) Beta-Bernoulli with Poisson(lambda) mixture - EB version

eb_PoissonBB_rare <- rarefaction(object = eb_PoissonBB_fit) # parameters of the Poisson's

# 4.2) Beta-Bernoulli with NB(n0,mu0) mixture - EB version

eb_NegBinBB_rare <- rarefaction(object = eb_NegBinBB_fit) # parameters of the NB's

# 4.3) IBP with Gamma(a,b) mixture - EB version

eb_GammaIBP_rare <- rarefaction(object = eb_GammaIBP_fit) # parameters of the NB's



# 5) Extrapolation ----

# 5.1) Beta-Bernoulli with Poisson(lambda) mixture - EB version
n <- nrow(feature_matrix_PoissonBB)
eb_PoissonBB_extr <- tibble(lambda = unname(unlist( 
  extrapolation(object = eb_PoissonBB_fit, M = 20)$lambda_post)),
  Kn = rep(Kn_PoissonBB, each = 20)) %>%
  mutate(lb = qpois(0.025, lambda, lower.tail = TRUE, log.p = FALSE),
         ub = qpois(0.975, lambda, lower.tail = TRUE, log.p = FALSE)) %>%
  rename(means = lambda) %>%
  add_row(means = 0, lb = 0, ub = 0, Kn = Kn_PoissonBB) %>%
  mutate(means = means + Kn, lb = lb + Kn, ub = ub + Kn) %>%
  add_column(x = c((n+1):(n+20), n),
             Model = "Poisson BB") %>%
  select(means, lb, ub, x, Model)


# 5.2) Beta-Bernoulli with NB(n0,mu0) mixture - EB version
n <- nrow(feature_matrix_NegBinBB)
eb_NegBinBB_extr <- tibble(mu0 = unname(unlist(
  extrapolation(object = eb_NegBinBB_fit, M = 20)$mu0_post )),
  n0 = unname(unlist( extrapolation(object = eb_NegBinBB_fit, M = 20)$n0_post )),
  Kn = rep(Kn_NegBinBB, each = 20)) %>%
  mutate(p = 1/(mu0/n0 + 1),
         lb = qnbinom(0.025, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE),
         ub = qnbinom(0.975, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE)) %>%
  rename(means = mu0) %>%
  add_row(means = 0, lb = 0, ub = 0, Kn = Kn_NegBinBB) %>%
  mutate(means = means + Kn, lb = lb + Kn, ub = ub + Kn) %>%
  add_column(x = c((n+1):(20+n), n),
             Model = paste0("NegBinomial BB x", eb_NegBinBB_fit$var_fct)) %>%
  select(means, lb, ub, x, Model)



# 5.3) IBP with Gamma(a,b) mixture - EB version
n <- nrow(feature_matrix_NegBinBB)
Kn_GammaIBP <- ncol(feature_matrix_GammaIBP[, colSums(feature_matrix_GammaIBP)!=0])
eb_GammaIBP_extr <- tibble(mu0 = unname(unlist(
  extrapolation(object = eb_GammaIBP_fit, M = 20)$mu0_post )),
  n0 = unname(unlist( extrapolation(object = eb_GammaIBP_fit, M = 20)$n0_post )),
  Kn = rep(Kn_GammaIBP, each = 20)) %>%
  mutate(p = 1/(mu0/n0 + 1),
         lb = qnbinom(0.025, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE),
         ub = qnbinom(0.975, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE)) %>%
  rename(means = mu0) %>%
  add_row(means = 0, lb = 0, ub = 0, Kn = Kn_GammaIBP) %>%
  mutate(means = means + Kn, lb = lb + Kn, ub = ub + Kn) %>%
  add_column(x = c((n+1):(20+n), n),
             Model = paste0("Gamma IBP, Variance: ", eb_GammaIBP_fit$var)) %>%
  select(means, lb, ub, x, Model)



