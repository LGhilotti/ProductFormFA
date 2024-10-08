#
# Application to Barro Colorado Island data ####
#

rm(list=ls())
library(ProductFormFA)
library(tidyverse)
library(scales)
library(ggpubr)
library(ggthemes)

source("R_script_paper/Routine_Chao.R")
source("R_script_paper/utils.R")

# Load data
library(vegan)
data("BCI")
data <- 1*as.matrix(BCI > 0)
data <- data[, colSums(is.na(data))==0]
data <- data[, colSums(data)!=0]

# Number of sites and number of species
n <- nrow(data)
Kn <- ncol(data)
print(paste0("Number of sites: ", n ))
print(paste0("Number of species: ", Kn))

# Randomly reorder sites
seed <- 12345
set.seed(seed)
data_mat <- data[sample.int(n, size = n, replace = F),]

# Plot accumulation
accum_df <- tibble( x = 0:n,
                    n_feat = c(0,rarefaction(data_mat[1:n,], n_reorderings = 1)))

ggplot(accum_df, aes(x = x, y = n_feat)) +
  geom_point(color="black", shape = 19, size = 0.1) + 
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) + 
  scale_color_tableau()
#ggsave(filename = "R_script_paper/Paper_plots/accumulation_BCI.pdf", width = 4, height = 4, dpi = 300, units = "in", device='pdf')




# EFPF approach -----

# Choices of variances
vars_fct_NegBinBB <- c(10, 1000)
vars_GammaIBP <- c(1, 1000)

# Initial parameters for optimization
eb_init_BB <- list(alpha = -10, s = 100, Nhat_prime = 200)
eb_known_BB <- list()

eb_init_IBP <- list(alpha = 0.5, s = 10, Gamma = 10)
eb_known_IBP <- list()

eb_params_obj_BB <- eb_params(model = "BB", 
                              init = eb_init_BB, known = eb_known_BB )
eb_params_obj_IBP <- eb_params(model = "IBP", 
                               init = eb_init_IBP, known = eb_known_IBP )

# PoissonBB
eb_EFPF_fit_PoissonBB <- GibbsFA_eb(feature_matrix = data_mat, 
                                    model = "PoissonBB", 
                                    type = "EFPF",
                                    eb_params =  eb_params_obj_BB)

# NegBinBB
list_eb_EFPF_fit_NegBinBB <-  vector(mode = "list", length = length(vars_fct_NegBinBB))
names(list_eb_EFPF_fit_NegBinBB) <- paste0("var_fct.", vars_fct_NegBinBB)

for (var_fct_NegBinBB in vars_fct_NegBinBB){
  
  list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]] <- 
    GibbsFA_eb(feature_matrix = data_mat,
               model = "NegBinBB", type = "EFPF",
               eb_params =  eb_params_obj_BB, 
               var_fct = var_fct_NegBinBB)
  
}

# GammaIBP
list_eb_EFPF_fit_GammaIBP <-  vector(mode = "list", length = length(vars_GammaIBP))
names(list_eb_EFPF_fit_GammaIBP) <- paste0("var.", vars_GammaIBP)

for (var_GammaIBP in vars_GammaIBP){
  
  list_eb_EFPF_fit_GammaIBP[[paste0("var.", var_GammaIBP)]] <-
    GibbsFA_eb(feature_matrix = data_mat,
               model = "GammaIBP", type = "EFPF",
               eb_params =  eb_params_obj_IBP,
               var_GammaIBP = var_GammaIBP)
  
}



## Model-checking on rarefaction ----
n_rare <- n
eb_EFPF_fit_PoissonBB_rare <- eb_EFPF_fit_PoissonBB
eb_EFPF_fit_NegBinBB_rare <- list_eb_EFPF_fit_NegBinBB[[1]]
eb_EFPF_fit_GammaIBP_rare <- list_eb_EFPF_fit_GammaIBP[[1]]

accum_df <- tibble( x = 0:n_rare,
                    n_feat = c(0,rarefaction(data_mat[1:n_rare,], n_reorderings = 50)))

rare_EFPF_PoissonBB <- tibble( lambda_post = unname(unlist(
  rarefaction(object = eb_EFPF_fit_PoissonBB_rare, seed = seed)$lambda_post ))) %>%
  rename(means = lambda_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "PoissonBB",
             x = c(1:n_rare,0))


rare_EFPF_NegBinBB <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_EFPF_fit_NegBinBB_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "NegBinBB",
             x = c(1:n_rare,0))



rare_EFPF_GammaIBP <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_EFPF_fit_GammaIBP_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "GammaIBP",
             x = c(1:n_rare,0))

# define quantities for plot
rare_EFPF_mixtureBB <- rare_EFPF_NegBinBB %>%
  mutate(Model = "PoissonBB/NegBinBB")

df_rare <- rbind(rare_EFPF_mixtureBB,
                 rare_EFPF_GammaIBP)

df_rare$Model <- factor(df_rare$Model,
                        levels = c("PoissonBB/NegBinBB", "GammaIBP"))

# for plot
ggplot(accum_df, aes(x = x, y = n_feat)) +
  geom_point(color="black", shape = 19, size = 0.8) + 
  geom_line(data = df_rare, aes(x = x, y = means, color = Model), linetype = "dashed") + 
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) + 
  scale_color_tableau(
    labels = c(
      "PoissonBB/NegBinBB" = "Mixtures of BBs",
      "GammaIBP" = "Mixtures of IBPs"
    )
  ) 
#ggsave(filename = "R_script_paper/Paper_plots/rarefaction_BCI_eb_EFPF.pdf", width = 4, height = 4, dpi = 300, units = "in", device='pdf')


## Model-checking on K_n_r -------
n_knr <- n
eb_EFPF_fit_PoissonBB_knr <- eb_EFPF_fit_PoissonBB
eb_EFPF_fit_NegBinBB_knr <- list_eb_EFPF_fit_NegBinBB[[1]]
eb_EFPF_fit_GammaIBP_knr <- list_eb_EFPF_fit_GammaIBP[[1]]

observed_K_n_r <- tibble( r = 1:n_knr,
                          k_n_r = K_n_r(data_mat[1:n_knr,], n_reorderings = 1)[[paste0('N = ', n_knr)]])


K_n_r_EFPF_PoissonBB <- tibble( lambda_est = unname(unlist(
  K_n_r(object = eb_EFPF_fit_PoissonBB_knr, n = n_knr)[[paste0('N = ', n_knr)]]$lambda_est ))) %>%
  rename(means = lambda_est) %>%
  add_column(Model = "PoissonBB",
             r = 1:n_knr)


K_n_r_EFPF_NegBinBB <- tibble( mu0_est = unname(unlist(
  K_n_r(object = eb_EFPF_fit_NegBinBB_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
  rename(means = mu0_est) %>%
  add_column(Model = "NegBinBB",
             r = 1:n_knr)


K_n_r_EFPF_GammaIBP <- tibble( mu0_est = unname(unlist(
  K_n_r(object = eb_EFPF_fit_GammaIBP_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
  rename(means = mu0_est) %>%
  add_column(Model = "GammaIBP",
             r = 1:n_knr)

# define quantities for plot
K_n_r_EFPF_mixtureBB <- K_n_r_EFPF_PoissonBB %>%
  mutate(Model = "PoissonBB/NegBinBB")

df_K_n_r <- rbind(K_n_r_EFPF_mixtureBB,
                  K_n_r_EFPF_GammaIBP)

df_K_n_r$Model <- factor(df_K_n_r$Model, 
                         levels = c("PoissonBB/NegBinBB", "GammaIBP"))

r_positive <- observed_K_n_r %>%
  filter(k_n_r > 0) %>%
  select(r) %>%
  filter(r < 15)

df_K_n_r_plot <- df_K_n_r %>%
  filter(r %in% c(r_positive$r))

observed_K_n_r_plot <- observed_K_n_r %>%
  filter(r %in% c(r_positive$r))

# for plot
ggplot(observed_K_n_r_plot,  aes(x = r, y = k_n_r)) +
  geom_point(color="black", shape = 19, size = 1) +
  geom_line( data = df_K_n_r_plot, aes(x = r, y = means, color = Model), linetype = "dashed") +
  scale_y_log10() +
  xlab("r") + ylab(expression(m[r])) + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) + 
  scale_color_tableau(
    labels = c(
      "PoissonBB/NegBinBB" = "Mixtures of BBs",
      "GammaIBP" = "Mixtures of IBPs"
    )
  )
#ggsave(filename = "R_script_paper/Paper_plots/knr_BCI_eb_EFPF.pdf", width = 4, height = 4, dpi = 300, units = "in", device='pdf')



## Prediction: richness and extrapolation -----------

### Richness -------

# PoissonBB
params_richness_EFPF_PoissonBB <- tibble( lambda_prime = 
                                            total_richness(eb_EFPF_fit_PoissonBB)$lambda_post) %>%
  add_column(Model = "Poisson BB") %>%
  mutate(lb = qpois(0.025, lambda_prime, lower.tail = TRUE, log.p = FALSE),
         ub = qpois(0.975, lambda_prime, lower.tail = TRUE, log.p = FALSE) )

# NegBinBB
params_richness_EFPF_NegBinBB <- tibble( n0_prime = numeric(),
                                         mu0_prime = numeric(),
                                         Model = character(),
                                         p_prime = numeric(),
                                         lb = numeric(), ub = numeric())

for (var_fct_NegBinBB in vars_fct_NegBinBB){
  
  eb_EFPF_NegBinBB_var <- list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]]
  
  params_richness_EFPF_NegBinBB_var <- tibble( n0_prime = 
                                                 total_richness(eb_EFPF_NegBinBB_var)$n0_post,
                                               mu0_prime = 
                                                 total_richness(eb_EFPF_NegBinBB_var)$mu0_post) %>%
    add_column(Model = paste0("NegBinomial BB x", var_fct_NegBinBB)) %>%
    mutate(p_prime = 1/(mu0_prime/n0_prime + 1),
           lb = qnbinom(0.025, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE),
           ub = qnbinom(0.975, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE) )
  
  
  params_richness_EFPF_NegBinBB <- bind_rows(params_richness_EFPF_NegBinBB,
                                             params_richness_EFPF_NegBinBB_var)
  
}



bounds <- tibble( lb = min(params_richness_EFPF_PoissonBB$lb + Kn, 
                           params_richness_EFPF_NegBinBB$lb + Kn),
                  ub = max(params_richness_EFPF_PoissonBB$ub + Kn,
                           params_richness_EFPF_NegBinBB$ub + Kn))

# PoissonBB
dens_richness_PoissonBB <- tibble( x = bounds$lb: bounds$ub) %>%
  mutate( y = dpois(x - Kn , lambda = params_richness_EFPF_PoissonBB$lambda_prime)) %>%
  add_column(Model = "Poisson BB")

# NegBinBB
dens_richness_NegBinBB <- tibble(x = integer(), y = numeric(),
                                 Model = character() )

for (var_fct_NegBinBB in vars_fct_NegBinBB){
  
  params_richness_EFPF_NegBinBB_var <- params_richness_EFPF_NegBinBB %>%
    filter(Model == paste0("NegBinomial BB x", var_fct_NegBinBB))
  
  dens_richness_NegBinBB_var <- tibble( x = bounds$lb : bounds$ub) %>%
    mutate( y = dnbinom(x - Kn, size = params_richness_EFPF_NegBinBB_var$n0_prime, 
                        prob = params_richness_EFPF_NegBinBB_var$p_prime)) %>%
    add_column(Model = paste0("NegBinomial BB x", var_fct_NegBinBB))
  
  dens_richness_NegBinBB <- bind_rows(dens_richness_NegBinBB,
                                      dens_richness_NegBinBB_var)
  
}



dens_richnesses <- rbind(dens_richness_PoissonBB,
                         dens_richness_NegBinBB) 


dens_richnesses$Model <- factor(dens_richnesses$Model, 
                                levels = c("Poisson BB",
                                           paste0("NegBinomial BB x", vars_fct_NegBinBB)))


# for plot
ggplot(dens_richnesses, aes(x = x, y = y, color = Model)) +
  geom_line() +
  theme_light() +
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  xlab("# distinct features") + rremove("ylab") +
  scale_color_tableau() +
  theme(aspect.ratio = 1)
#ggsave(filename = "R_script_paper/Paper_plots/richness_BCI_eb_EFPF.pdf", width = 5, height = 5, dpi = 300, units = "in", device='pdf')


# Compute mean and variance species richness

# PoissonBB
rich_pars <- params_richness_EFPF_PoissonBB 

print(paste0("mean richness = ", rich_pars$lambda_prime + Kn ))
Kn + qpois(c(0.025, 0.975), lambda = rich_pars$lambda_prime)


# NegBinBB
rich_pars <- params_richness_EFPF_NegBinBB %>%
  filter(Model == "NegBinomial BB x10")

print(paste0("mean richness = ", rich_pars$mu0_prime + Kn ))
Kn + qnbinom(c(0.025, 0.975), size = rich_pars$n0_prime, prob = rich_pars$p_prime)


### Extrapolation -----

# Extract accumulation curve of the observed sample (or average accumulation)
M <- 400

accum_df <- tibble( x = 0:n,
                    n_feat = c(0,rarefaction(data_mat, n_reorderings = 200)))


# PoissonBB
extr_EFPF_PoissonBB_df <- tibble(lambda = unname(unlist( 
  extrapolation(object = eb_EFPF_fit_PoissonBB, M = M, seed = seed)$lambda_post)),
  Kn = rep(Kn, each = M)) %>%
  mutate(lb = qpois(0.025, lambda, lower.tail = TRUE, log.p = FALSE),
         ub = qpois(0.975, lambda, lower.tail = TRUE, log.p = FALSE)) %>%
  rename(means = lambda) %>%
  add_row(means = 0, lb = 0, ub = 0, Kn = Kn) %>%
  mutate(means = means + Kn, lb = lb + Kn, ub = ub + Kn) %>%
  add_column(x = c((n+1):(n+M), n),
             Model = "Poisson BB")

extr_EFPF_PoissonBB_df$x <- as.integer(extr_EFPF_PoissonBB_df$x)
extr_EFPF_PoissonBB_df <- extr_EFPF_PoissonBB_df %>%
  select(means, lb, ub, x, Model)

# NegBin
extr_EFPF_NegBinBB_df <- tibble(means = numeric(), 
                                lb = numeric(), ub = numeric(),
                                x = integer(), Model = character())

for (var_fct_NegBinBB in vars_fct_NegBinBB){
  
  eb_EFPF_NegBinBB_var <- list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]]
  
  extr_EFPF_NegBinBB_df_var <- tibble(mu0 = unname(unlist(
    extrapolation(object = eb_EFPF_NegBinBB_var, M = M, seed = seed)$mu0_post )),
    n0 = unname(unlist( extrapolation(object = eb_EFPF_NegBinBB_var, M = M, seed = seed)$n0_post )),
    Kn = rep(Kn, each = M)) %>%
    mutate(p = 1/(mu0/n0 + 1),
           lb = qnbinom(0.025, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE),
           ub = qnbinom(0.975, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE)) %>%
    rename(means = mu0) %>%
    add_row(means = 0, lb = 0, ub = 0, Kn = Kn) %>%
    mutate(means = means + Kn, lb = lb + Kn, ub = ub + Kn) %>%
    add_column(x = c((n+1):(M+n), n),
               Model = paste0("NegBinomial BB x", var_fct_NegBinBB))
  
  extr_EFPF_NegBinBB_df_var$x <- as.integer(extr_EFPF_NegBinBB_df_var$x)
  extr_EFPF_NegBinBB_df_var <- extr_EFPF_NegBinBB_df_var %>%
    select(means, lb, ub, x, Model)
  
  extr_EFPF_NegBinBB_df <- bind_rows(extr_EFPF_NegBinBB_df, 
                                     extr_EFPF_NegBinBB_df_var)
}


extr_EFPF_PoissonBB_df <- extr_EFPF_PoissonBB_df %>%
  add_column(Model_gen = "PoissonBB/NegBinBB")
extr_EFPF_NegBinBB_df <- extr_EFPF_NegBinBB_df %>%
  add_column(Model_gen = "PoissonBB/NegBinBB") 

extr_all_df <- rbind(extr_EFPF_PoissonBB_df, 
                     extr_EFPF_NegBinBB_df)

extr_all_df$Model <- factor(extr_all_df$Model,
                            levels = c("Poisson BB",
                                       paste0("NegBinomial BB x", vars_fct_NegBinBB)))

extr_all_df$Model_gen <- factor(extr_all_df$Model_gen,
                                levels = c("PoissonBB/NegBinBB", "GammaIBP"))

# for plot
ggplot(extr_all_df, aes(x, means, color = Model)) +
  geom_line(linetype = "dashed") +
  geom_point( data = accum_df, aes(x = x, y = n_feat),
              color="black", shape = 19, size = 0.1) +
  geom_ribbon(aes(ymin = lb, ymax = ub, color = Model), alpha = 0) +
  geom_vline(aes(xintercept = n) , linetype = "dashed", color = "grey") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) +
  scale_color_tableau()
#ggsave(filename = "R_script_paper/Paper_plots/extr_BCI_eb_EFPF.pdf", width = 5.2, height = 5.2, dpi = 300, units = "in", device='pdf')


# Compute extrapolation on a grid: numerical values
extr_EFPF_NegBinBB_df %>%
  filter(Model == "NegBinomial BB x10",
         x %in% c(n + 1, n + 10, n + 100, n + 1000)) %>%
  mutate(means_new = means - Kn,
         lb_new = lb -Kn,
         ub_new = ub - Kn)

extr_EFPF_PoissonBB_df %>%
  filter(x %in% c(n + 1, n + 10, n + 100, n + 1000)) %>%
  mutate(means_new = means - Kn,
         lb_new = lb -Kn,
         ub_new = ub - Kn)



# Fully-Bayesian approach ----------

# We focus on NegBinBB's + priors (since it is selected from model-checking)

# Fit and estimate richness, rarefaction and extrapolation for GibbsFA's (save workspace)
if (!file.exists("R_script_paper/fullybayes_BCI_fit_estimate.RData")) {
  
  list_prior_fit_NegBinBB <-  vector(mode = "list", length = length(vars_fct_NegBinBB))
  names(list_prior_fit_NegBinBB) <- paste0("var_fct.", vars_fct_NegBinBB)
  
  # Initialization and MCMC setting 
  mcmcparams_NegBinBB <- list(tau = 0.1, 
                              S = 5*10^4, n_burnin = 5*10^3, thin = 2)
  mcmcparams_obj_NegBinBB <- mcmcparameters(model = "NegBinBB", mcmcparams = mcmcparams_NegBinBB)
  
  init_NegBinBB <- list(alpha_0 = - 1, s_0 = 15)
  init_obj_NegBinBB <- initialization(model = "NegBinBB", init = init_NegBinBB )
  
  # EB estimates
  small_val <- 10^(-3)
  alpha_eb <- list_eb_EFPF_fit_NegBinBB[[1]]$alpha
  theta_eb <- list_eb_EFPF_fit_NegBinBB[[1]]$theta
  
  s_eb <- alpha_eb + theta_eb
  
  print(paste0("Prior variance of -alpha: ", - alpha_eb/small_val  ))
  print(paste0("Prior variance of s: ", s_eb/small_val  ))
  
  # Fit the model
  for (var_fct in vars_fct_NegBinBB){
    
    n0_eb <- list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct)]]$n0
    mu0_eb <- list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct)]]$mu0
    
    hyper_NegBinBB <- list(a_alpha = - alpha_eb*small_val, b_alpha = small_val,
                           a_s = s_eb*small_val , b_s = small_val,
                           n0 = n0_eb, mu0 = mu0_eb)
    prior_obj_NegBinBB <- prior(model = "NegBinBB", hyper = hyper_NegBinBB)
    
    
    list_prior_fit_NegBinBB[[paste0("var_fct.", var_fct)]] <- 
      GibbsFA(feature_matrix = data_mat,
              model = "NegBinBB", 
              prior = prior_obj_NegBinBB,
              initialization = init_obj_NegBinBB,
              mcmcparams = mcmcparams_obj_NegBinBB)
    
  }
  
  # Save the entire workspace related to the type just performed
  save(list = ls(all.names = TRUE), file =  "R_script_paper/fullybayes_BCI_fit_estimate.RData")
  
}

# Load the Work space
load("R_script_paper/fullybayes_BCI_fit_estimate.RData")


## Convergence checks --------
library(ggmcmc)
library(coda)

params_prior_NegBinBB <- list_prior_fit_NegBinBB[[paste0("var_fct.", vars_fct_NegBinBB[1])]][c("n0_chain","mu0_chain", "alpha_chain", "theta_chain")]
params_prior_NegBinBB_df <- as.data.frame(do.call(cbind, params_prior_NegBinBB))

samples_NegBinBB <- mcmc.list(mcmc(params_prior_NegBinBB_df))
samples_ggs_NegBinBB <- ggs(samples_NegBinBB, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_NegBinBB)

effectiveSize(params_prior_NegBinBB_df)


## Prediction: richness and extrapolation -----------

### Extrapolation ------

if (!file.exists("R_script_paper/fullybayes_BCI_extrapolation.RData")) {
  
  M <- 400
  
  # NegBinBB
  extr_prior_NegBinBB_df <- tibble(means = numeric(), 
                                   lb = numeric(), ub = numeric(),
                                   x = integer(), Model = character())
  
  
  
  for (var_fct in vars_fct_NegBinBB){
    
    extr_prior_NegBinBB_var <- extrapolation(object = list_prior_fit_NegBinBB[[paste0("var_fct.", var_fct)]],
                                             M = M) 
    
    extr_prior_NegBinBB_var_df_tmp <- as_tibble(t(bind_rows(as.data.frame(lapply(extr_prior_NegBinBB_var, quantile, prob = c(0.025, 0.975))),
                                                            as.data.frame(lapply(extr_prior_NegBinBB_var, mean))))) 
    colnames(extr_prior_NegBinBB_var_df_tmp) <- c("lb", "ub", "means")
    
    extr_prior_NegBinBB_var_df <- extr_prior_NegBinBB_var_df_tmp %>%
      add_column(x = 1:nrow(extr_prior_NegBinBB_var_df_tmp),
                 Model = paste0("NegBinBB x", var_fct)) %>%
      mutate( x = x + n ) %>%
      add_row(means = Kn, ub = Kn, lb = Kn, x=n, 
              Model = paste0("NegBinBB x", var_fct))
    
    extr_prior_NegBinBB_var_df$x <- as.integer(extr_prior_NegBinBB_var_df$x)
    
    extr_prior_NegBinBB_var_df <- extr_prior_NegBinBB_var_df %>%
      select(means, lb, ub, x, Model)
    
    extr_prior_NegBinBB_df <- bind_rows(extr_prior_NegBinBB_df, 
                                        extr_prior_NegBinBB_var_df)
  }
  
  
  # Save the entire workspace related to the type just performed
  save(list = ls(all.names = TRUE), file =  "R_script_paper/fullybayes_BCI_extrapolation.RData")
  
}

# Load the Work space
load("R_script_paper/fullybayes_BCI_extrapolation.RData")
extr_prior_NegBinBB_df_final <- extr_prior_NegBinBB_df %>%
  mutate(Model = case_when(
    Model == "NegBinBB x10" ~  "NegBinomial BB x10",
    Model == "NegBinBB x1000" ~ "NegBinomial BB x1000")) %>%
  add_column(Type = "Fully Bayesian")

# Join the df related to prior and EFPF to compare in the plot
extr_EFPF_NegBinBB_df_final <- extr_EFPF_NegBinBB_df %>%
  add_column(Type = "EB")


extr_joint_NegBinBB_df <- bind_rows(extr_EFPF_NegBinBB_df_final,
                                    extr_prior_NegBinBB_df_final)


# for plot
ggplot(extr_joint_NegBinBB_df, aes(x, means, color = Type )) +
  geom_line(linetype = "dashed") +
  facet_wrap(. ~ Model,  scales = "free_x") +
  geom_point( data = accum_df, aes(x = x, y = n_feat),
              color="black", shape = 19, size = 0.1) +
  geom_ribbon(aes(ymin = lb, ymax = ub, color = Type), alpha = 0) +
  geom_vline(aes(xintercept = n) , linetype = "dashed", color = "grey") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  labs(color = "Approach") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) +
  scale_color_tableau()
#ggsave(filename = "R_script_paper/Paper_plots/extr_BCI_fullybayes.pdf", width = 6, height = 4, dpi = 300, units = "in", device='pdf')




### Richness ------

bounds <- list("lb" = 250, "ub" = 400)

# for EFPF approach
richness_EFPF_NegBinBB_df <- tibble(x = integer(), y = numeric(),
                                    Model = character() )

for (var_fct_NegBinBB in vars_fct_NegBinBB){
  
  params_richness_EFPF_NegBinBB_var <- params_richness_EFPF_NegBinBB %>%
    filter(Model == paste0("NegBinomial BB x", var_fct_NegBinBB))
  
  dens_richness_NegBinBB_var <- tibble( x = bounds$lb : bounds$ub) %>%
    mutate( y = dnbinom(x - Kn, size = params_richness_EFPF_NegBinBB_var$n0_prime, 
                        prob = params_richness_EFPF_NegBinBB_var$p_prime)) %>%
    add_column(Model = paste0("NegBinomial BB x", var_fct_NegBinBB))
  
  richness_EFPF_NegBinBB_df <- bind_rows(richness_EFPF_NegBinBB_df,
                                         dens_richness_NegBinBB_var)
  
}

richness_EFPF_NegBinBB_df$Model <- factor(richness_EFPF_NegBinBB_df$Model, 
                                           levels = paste0("NegBinomial BB x", vars_fct_NegBinBB))


# for Fully-Bayes approach
richness_prior_NegBinBB_df <- tibble(Model = character(),
                                     y = numeric())
           
for (var_fct_NegBinBB in vars_fct_NegBinBB){

  prior_NegBinBB_var <- list_prior_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]]
  
  richness_prior_NegBinBB_df_var <- tibble(
    y = total_richness(object = prior_NegBinBB_var),
    Model = paste0("NegBinomial BB x", var_fct_NegBinBB)
  ) 
  
  richness_prior_NegBinBB_df <- bind_rows(richness_prior_NegBinBB_df,
                                          richness_prior_NegBinBB_df_var)

}


richness_prior_NegBinBB_df$Model <- factor(richness_prior_NegBinBB_df$Model, 
                                           levels = paste0("NegBinomial BB x", vars_fct_NegBinBB))

# prepare final df
richness_EFPF_NegBinBB_df_final <- richness_EFPF_NegBinBB_df %>%
  add_column(Type = "EB")

richness_prior_NegBinBB_df_final <- richness_prior_NegBinBB_df %>%
  add_column(Type = "Fully Bayesian")


# for plot
ggplot() +
  geom_line(data = richness_EFPF_NegBinBB_df_final, aes(x = x, y = y, color = Type) ) +
  stat_density(data = richness_prior_NegBinBB_df_final, aes(x=y, color = Type), geom="line",position="identity", bw = 3) +
  theme_light() +
  facet_wrap(~Model) +
  theme(legend.position = "top") +
  labs(color = "Approach") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(
    limits = c(250, 400) 
    ) +
  xlab("# distinct features") + rremove("ylab") +
  scale_color_tableau() +
  theme(aspect.ratio = 1)
#ggsave(filename = "R_script_paper/Paper_plots/richness_BCI_prior.pdf", width = 6, height = 4, dpi = 300, units = "in", device='pdf')


# Compute mean and variance of richness

rich_draws <- richness_prior_NegBinBB_df %>%
  filter(Model == "NegBinomial BB x1000")

print(paste0("mean of N: ", mean(rich_draws$y)))
quantile(rich_draws$y, prob = c(0.025, 0.975))

