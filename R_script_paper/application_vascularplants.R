#
# Application to Vascular plants data ####
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
data <- read.csv(file = "R_script_paper/mazz2016_Plants.csv", header = TRUE,
                 row.names="X")
data <- t(data)
data <- data[, colSums(is.na(data))==0]
data <- data[, colSums(data)!=0]

# Number of sites and number of species
n <- nrow(data)
Kn <- ncol(data)
print(paste0("Number of sites: ", n ))
print(paste0("Number of species: ", Kn))

# Randomly reorder sites
seed <- 12345678
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
#ggsave(filename = "R_script_paper/Paper_plots/accumulation_Plants.pdf", width = 4, height = 4, dpi = 300, units = "in", device='pdf')



# EFPF approach -----

# Choices of variances
vars_fct_NegBinBB <- c(10, 1000)
vars_GammaIBP <- c(0.01, 100)

# Initial parameters for optimization
eb_init_BB <- list(alpha = -10, s = 100, Nhat_prime = 200)
eb_known_BB <- list()

eb_init_IBP <- list(alpha = 0.4, s = 2, Gamma = 10)
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
                    n_feat = c(0,rarefaction(data_mat[1:n_rare,], n_reorderings = 200)))

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
  geom_point(color="black", shape = 18, size = 0.8) + 
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
#ggsave(filename = "R_script_paper/Paper_plots/rarefaction_Plants_eb_EFPF.pdf", width = 4, height = 4, dpi = 300, units = "in", device='pdf')


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
  filter(r < 9)

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
#ggsave(filename = "R_script_paper/Paper_plots/knr_Plants_eb_EFPF.pdf", width = 4, height = 4, dpi = 300, units = "in", device='pdf')



## Extrapolation -----------

# 0) alpha-diversity of the gamma mixture
idx_gamma <- 2 # or 1, for the 2 different choices of the prior variance of gamma

a <- list_eb_EFPF_fit_GammaIBP[[idx_gamma]]$a
b <- list_eb_EFPF_fit_GammaIBP[[idx_gamma]]$b
alpha <- list_eb_EFPF_fit_GammaIBP[[idx_gamma]]$alpha
theta <- list_eb_EFPF_fit_GammaIBP[[idx_gamma]]$theta

gamma_a_t_n <- sum(exp(lgamma(alpha + theta + (1:n) - 1) - lgamma(alpha + theta) -
                         lgamma(theta + (1:n)) + lgamma(theta +1) ) )
# posterior diversity
a_gamma <- a + Kn
b_gamma <- (b + gamma_a_t_n)*gamma(theta+alpha)/gamma(theta+1)*alpha
print(paste0("mean(diversity) = ", a_gamma/b_gamma))  
print(paste0("var(diversity) = ", a_gamma/b_gamma^2))  

# Extract accumulation curve of the observed sample (or average accumulation)

M <- 400

accum_df <- tibble( x = 0:n,
                    n_feat = c(0,rarefaction(data_mat, n_reorderings = 200)))


# GammaIBP
extr_EFPF_GammaIBP_df <- tibble(means = numeric(), 
                                lb = numeric(), ub = numeric(),
                                x = integer(), Model = character())

for (var_GammaIBP in vars_GammaIBP){
  
  eb_EFPF_GammaIBP_var <- list_eb_EFPF_fit_GammaIBP[[paste0("var.", var_GammaIBP)]]
  
  extr_GammaIBP_df_var <- tibble(mu0 = unname(unlist(
    extrapolation(object = eb_EFPF_GammaIBP_var, M = M, seed = seed)$mu0_post )),
    n0 = unname(unlist( extrapolation(object = eb_EFPF_GammaIBP_var, M = M, seed = seed)$n0_post )),
    Kn = rep(Kn, each = M)) %>%
    mutate(p = 1/(mu0/n0 + 1),
           lb = qnbinom(0.025, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE),
           ub = qnbinom(0.975, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE)) %>%
    rename(means = mu0) %>%
    add_row(means = 0, lb = 0, ub = 0, Kn = Kn) %>%
    mutate(means = means + Kn, lb = lb + Kn, ub = ub + Kn) %>%
    add_column(x = c((n+1):(M+n), n),
               Model = paste0("Gamma IBP, Variance: ", var_GammaIBP))
  
  extr_GammaIBP_df_var$x <- as.integer(extr_GammaIBP_df_var$x)
  extr_GammaIBP_df_var <- extr_GammaIBP_df_var %>%
    select(means, lb, ub, x, Model)
  
  extr_EFPF_GammaIBP_df <- bind_rows(extr_EFPF_GammaIBP_df, 
                                     extr_GammaIBP_df_var)
}

extr_EFPF_GammaIBP_df <- extr_EFPF_GammaIBP_df %>%
  add_column(Model_gen = "GammaIBP")
extr_all_df <- rbind(extr_EFPF_GammaIBP_df)


extr_all_df$Model <- factor(extr_all_df$Model,
                            levels = c("Poisson BB",
                                       paste0("NegBinomial BB x", vars_fct_NegBinBB),
                                       paste0("Gamma IBP, Variance: ", vars_GammaIBP)))

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
#ggsave(filename = "R_script_paper/Paper_plots/extr_Plants_eb_EFPF.pdf", width = 5, height = 5, dpi = 300, units = "in", device='pdf')



# Compute extrapolation on a grid: numerical values

extr_EFPF_GammaIBP_df %>%
  filter(Model == "Gamma IBP, Variance: 0.01",
         x %in% c(n + 1, n + 10, n + 100, n + 1000)) %>%
  mutate(means_new = means - Kn,
         lb_new = lb -Kn,
         ub_new = ub - Kn)



# Fully-Bayesian approach ----------

# We focus on GammaIBP + prior (since it is selected from model-checking)
  
# Fit for GibbsFA's (save workspace)
if (!file.exists("R_script_paper/fullybayes_Plants_fit_estimate.RData")) {
  
  list_prior_fit_GammaIBP <-  vector(mode = "list", length = length(vars_GammaIBP))
  names(list_prior_fit_GammaIBP) <- paste0("var.", vars_GammaIBP)
  
  # Initialization and MCMC setting 
  mcmcparams_GammaIBP <- list(sigq_alpha = 0.1, sigq_s = 0.1, 
                              S = 5*10^4, n_burnin = 5*10^3, thin = 2)
  mcmcparams_obj_GammaIBP <- mcmcparameters(model = "GammaIBP", mcmcparams = mcmcparams_GammaIBP)
  
  init_GammaIBP <- list(alpha_0 = 0.5, s_0 = 15)
  init_obj_GammaIBP <- initialization(model = "GammaIBP_single_prior", init = init_GammaIBP )
  
  # EB estimates
  small_val_alpha <- 2
  small_val_s <- 10^(-2)
  
  alpha_eb <- list_eb_EFPF_fit_GammaIBP[[1]]$alpha
  theta_eb <- list_eb_EFPF_fit_GammaIBP[[1]]$theta
  
  t_eb <- (1 - alpha_eb)/alpha_eb
  s_eb <- alpha_eb + theta_eb
  
  print(paste0("Prior variance of alpha: ", 
               t_eb/(1 + t_eb)^2 /(1 + small_val_alpha*(1+t_eb))))
  
  print(paste0("Prior variance of s: ", 
               s_eb/small_val_s))
  
  # Fit the model
  for (var_GammaIBP in vars_GammaIBP){
    
    a_eb <- list_eb_EFPF_fit_GammaIBP[[paste0("var.", var_GammaIBP)]]$a
    b_eb <- list_eb_EFPF_fit_GammaIBP[[paste0("var.", var_GammaIBP)]]$b
    
    # Hyperparameters elicitation 
    hyper_GammaIBP <- list(a = a_eb, b = b_eb,
                           a_alpha = small_val_alpha, b_alpha = t_eb*small_val_alpha,
                           a_s = s_eb*small_val_s , b_s = small_val_s)
    prior_obj_GammaIBP <- prior(model = "GammaIBP_single_prior", hyper = hyper_GammaIBP)
    
    
    list_prior_fit_GammaIBP[[paste0("var.", var_GammaIBP)]] <- 
      GibbsFA(feature_matrix = data_mat,
              model = "GammaIBP_single_prior", 
              prior = prior_obj_GammaIBP,
              initialization = init_obj_GammaIBP,
              mcmcparams = mcmcparams_obj_GammaIBP)
    
  }
  
  # Save the entire workspace related to the type just performed
  save(list = ls(all.names = TRUE), file =  "R_script_paper/fullybayes_Plants_fit_estimate.RData")
  
}


# Load the Work space
load("R_script_paper/fullybayes_Plants_fit_estimate.RData")


## Convergence checks --------
library(ggmcmc)
library(coda)

params_prior_GammaIBP <- list_prior_fit_GammaIBP[[paste0("var.", vars_GammaIBP[2])]][c("a_chain","b_chain", "alpha_chain", "theta_chain")]
params_prior_GammaIBP_df <- as.data.frame(do.call(cbind, params_prior_GammaIBP))

samples_GammaIBP <- mcmc.list(mcmc(params_prior_GammaIBP_df))
samples_ggs_GammaIBP <- ggs(samples_GammaIBP, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_GammaIBP)

effectiveSize(params_prior_GammaIBP_df)


## Extrapolation ------

if (!file.exists("R_script_paper/fullybayes_Plants_extrapolation.RData")) {
  
  M <- 400
  
  # GammaIBP
  extr_prior_GammaIBP_df <- tibble(means = numeric(), 
                                  lb = numeric(), ub = numeric(),
                                  x = integer(), Model = character())
  
  for (var_GammaIBP in vars_GammaIBP){
    
    extr_prior_GammaIBP_var <- extrapolation(object = list_prior_fit_GammaIBP[[paste0("var.", var_GammaIBP)]],
                                             M = M) 
    
    extr_prior_GammaIBP_var_df_tmp <- as_tibble(t(bind_rows(as.data.frame(lapply(extr_prior_GammaIBP_var, quantile, prob = c(0.025, 0.975))),
                                                    as.data.frame(lapply(extr_prior_GammaIBP_var, mean))))) 
    colnames(extr_prior_GammaIBP_var_df_tmp) <- c("lb", "ub", "means")
    
    extr_prior_GammaIBP_var_df <- extr_prior_GammaIBP_var_df_tmp %>%
      add_column(x = 1:nrow(extr_prior_GammaIBP_var_df_tmp),
                 Model = paste0("GammaIBP, var = ", var_GammaIBP)) %>%
      mutate( x = x + n ) %>%
      add_row(means = Kn, ub = Kn, lb = Kn, x=n, 
              Model = paste0("GammaIBP, var = ", var_GammaIBP))
    
    extr_prior_GammaIBP_var_df$x <- as.integer(extr_prior_GammaIBP_var_df$x)
    
    extr_prior_GammaIBP_var_df <- extr_prior_GammaIBP_var_df %>%
      select(means, lb, ub, x, Model)
    
    extr_prior_GammaIBP_df <- bind_rows(extr_prior_GammaIBP_df, 
                                        extr_prior_GammaIBP_var_df)
  }
  
  # Save the entire workspace related to the type just performed
  save(list = ls(all.names = TRUE), file =  "R_script_paper/fullybayes_Plants_extrapolation.RData")
  
}

# Load the Work space
load("R_script_paper/fullybayes_Plants_extrapolation.RData")
extr_prior_GammaIBP_df_final <- extr_prior_GammaIBP_df %>%
  mutate(Model = case_when(
    Model == "GammaIBP, var = 0.01" ~  "Gamma IBP, Variance: 0.01",
    Model == "GammaIBP, var = 100" ~ "Gamma IBP, Variance: 100")) %>%
  add_column(Type = "Fully Bayesian")

# Join the df related to prior and EFPF to compare in the plot
extr_EFPF_GammaIBP_df_final <- extr_EFPF_GammaIBP_df %>%
  add_column(Type = "EB")


extr_joint_GammaIBP_df <- bind_rows(extr_EFPF_GammaIBP_df_final,
                                    extr_prior_GammaIBP_df_final)


# for plot
ggplot(extr_joint_GammaIBP_df, aes(x, means, color = Type )) +
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
#ggsave(filename = "R_script_paper/Paper_plots/extr_Plants_fullybayes.pdf", width = 6, height = 4, dpi = 300, units = "in", device='pdf')

