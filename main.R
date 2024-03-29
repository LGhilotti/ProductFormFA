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

sample_GammaIBP <- rGammaIBP(alpha = 0.5, theta = 10, a = 1, b = 1, n = 20)

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
mcmcparams <- list(tau = 0.1, S = 200, n_burnin = 10, thin = 2)

mcmcparams_obj <- mcmcparameters(model = "PoissonBB", mcmcparams = mcmcparams)

# Fit the model
PoissonBB_fit <- GibbsFA(feature_matrix = feature_matrix, 
                       model = "PoissonBB", 
                       prior = prior_obj, 
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
                        prior = prior_obj, 
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

# Fit the model
GammaIBP_fit <- GibbsFA(feature_matrix = feature_matrix, 
                        model = "GammaIBP", 
                        prior = prior_obj, 
                        initialization = init_obj, 
                        mcmcparams = mcmcparams_obj)


# 3) Estimate the total richness ----

# 3.1) Beta-Bernoulli with Poisson(lambda) mixture

PoissonBB_rich <- total_richness(object = PoissonBB_fit)

# 3.2) Beta-Bernoulli with NB(n0,mu0) mixture

NegBinBB_rich <- total_richness(object = NegBinBB_fit)


# 4) Extrapolation ----

# 4.1) Beta-Bernoulli with Poisson(lambda) mixture

PoissonBB_extr <- extrapolation(object = PoissonBB_fit, M = 20) # option: only_last = TRUE, seed

# 4.2) Beta-Bernoulli with NB(n0,mu0) mixture

NegBinBB_extr <- extrapolation(object = NegBinBB_fit, M = 20) # option: only_last = TRUE, seed

# 4.3) IBP with Gamma(a,b) mixture

GammaIBP_extr <- extrapolation(object = GammaIBP_fit, M = 20) # option: only_last = TRUE, seed


# 5) Rarefaction ----

# 5.1) Beta-Bernoulli with Poisson(lambda) mixture

PoissonBB_rare <- rarefaction(object = PoissonBB_fit) # option: seed

# 5.2) Beta-Bernoulli with NB(n0,mu0) mixture

NegBinBB_rare <- rarefaction(object = NegBinBB_fit) # option: seed

# 5.3) IBP with Gamma(a,b) mixture

GammaIBP_rare <- rarefaction(object = GammaIBP_fit) # option: seed



# 6) Plot richness ----

# 6.1) Beta-Bernoulli with Poisson(lambda) mixture

plot(x = PoissonBB_fit, type = "richness")

# 6.2) Beta-Bernoulli with NB(n0,mu0) mixture

plot(x = NegBinBB_fit, type = "richness")


# 7) Plot Rarefaction ----

# 7.1) Beta-Bernoulli with Poisson(lambda) mixture

plot(x = PoissonBB_fit, type = "rarefaction") # option: plot_sample = FALSE , seed

# 7.2) Beta-Bernoulli with NB(n0,mu0) mixture

plot(x = NegBinBB_fit, type = "rarefaction") # option: plot_sample = FALSE , seed

# 7.3) IBP with Gamma(a,b) mixture

plot(x = GammaIBP_fit, type = "rarefaction") # option: plot_sample = FALSE , seed


# 8) Plot Extrapolation ----

# 8.1) Beta-Bernoulli with Poisson(lambda) mixture

plot(x = PoissonBB_fit, type = "extrapolation", M = 20 ) # option: seed

# 8.2) Beta-Bernoulli with NB(n0,mu0) mixture

plot(x = NegBinBB_fit, type = "extrapolation", M = 20) # option:  seed

# 8.3) IBP with Gamma(a,b) mixture

plot(x = GammaIBP_fit, type = "extrapolation", M = 20) # option:  seed
