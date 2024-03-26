#
#### How to use the functions? ------
#
rm(list = ls())
library(ProductFormFA)

# 1) Random generation of the processes ----

seed = 1234

# 1.1) Beta-Bernoulli with Poisson(lambda) mixture

sample_PoissonBB <- rPoissonBB(alpha = -1, theta = 2, lambda = 10, n = 10)

sample_PoissonBB

# 1.2) Beta-Bernoulli with NB(n0, mu0) mixture

sample_NegBinBB <- rNegBinBB(alpha = -1, theta = 2, n0 = 10, mu0 = 0.5, n = 10)

sample_NegBinBB

# 1.3) IBP with Gamma(a, b) mixture

sample_GammaIBP <- rGammaIBP(alpha = 0.5, theta = 2, a = 1, b = 1, n = 10)

sample_GammaIBP


# 2) Run the models on generated data ----

# 2.1) Beta-Bernoulli with Poisson(lambda) mixture

feature_matrix <- convert_features_matrix(sample_PoissonBB$features)

# Hyperparameters
hyper <- list(a_alpha = 1, b_alpha = 0.1,
              a_s = 2, b_s = 0.2,
              lambda = 10)

prior_obj <- prior(model = "PoissonBB", hyper = hyper) 
                  
summary(prior_obj)

# Initialization 
init <- list(alpha_0 = -1, s_0 = 1)

init_obj <- initialization(model = "PoissonBB", init = init )

# MCMC setting
mcmcparams <- list(tau = 0.1, S = 100, n_burnin = 10, thin = 2)

mcmcparams_obj <- mcmcparameters(model = "PoissonBB", mcmcparams = mcmcparams)

# Fit the model
PoissonBB_fit <- GibbsFA(feature_matrix = feature_matrix, 
                       model = "PoissonBB", 
                       hyperparams = prior_obj, 
                       initialization = init_obj, 
                       mcmcparams = mcmcparams_obj)


# 2.2) Beta-Bernoulli with NB(n0, mu0) mixture

feature_matrix <- convert_features_matrix(sample_NegBinBB$features)

# Hyperparameters
hyper <- list(a_alpha = 1, b_alpha = 0.1,
              a_s = 2, b_s = 0.2,
              n0 = 10, mu0 = 0.5)

prior_obj <- prior(model = "NegBinBB", hyper = hyper) 

summary(prior_obj)

# Initialization 
init <- list(alpha_0 = -1, s_0 = 1)

init_obj <- initialization(model = "NegBinBB", init = init )

# MCMC setting
mcmcparams <- list(tau = 0.1, S = 100, n_burnin = 10, thin = 2)

mcmcparams_obj <- mcmcparameters(model = "NegBinBB", mcmcparams = mcmcparams)

# Fit the model
NegBinBB_fit <- GibbsFA(feature_matrix = feature_matrix, 
                        model = "NegBinBB", 
                        hyperparams = prior_obj, 
                        initialization = init_obj, 
                        mcmcparams = mcmcparams_obj)


# 2.3) IBP with Gamma(a,b) mixture

feature_matrix <- convert_features_matrix(sample_GammaIBP$features)

# Hyperparameters
hyper <- list(a_alpha = 1, b_alpha = 2,
                    a_s = 2, b_s = 0.2,
                    q = 0.5, r = 1, t = 1)

prior_obj <- prior(model = "GammaIBP", hyper = hyper) 

summary(prior_obj)

# Initialization 
init <- list(alpha_0 = 0.5, s_0 = 1, a_0 = 1, b_0 = 1)

init_obj <- initialization(model = "GammaIBP", init = init )

# MCMC setting
mcmcparams <- list(sigq_alpha = 0.1, sigq_s = 1, 
                   S = 100, n_burnin = 10, thin = 2)

mcmcparams_obj <- mcmcparameters(model = "GammaIBP", mcmcparams = mcmcparams)

GammaIBP_fit <- GibbsFA(feature_matrix = feature_matrix, 
                        model = "GammaIBP", 
                        hyperparams = prior_obj, 
                        initialization = init_obj, 
                        mcmcparams = mcmcparams_obj)
