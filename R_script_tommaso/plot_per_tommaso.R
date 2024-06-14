rm(list=ls())
library(ProductFormFA)
library(tidyverse)
library(scales)
library(ggpubr)
library(ggthemes)
library(patchwork)
library(latex2exp)
library(dplyr, warn.conflicts = FALSE)

# Load the Work space
mechanism <- "custom"
load(paste0("R_script_tommaso/eb_EFPF_",mechanism,"_fit_estimate_singledataset.RData"))
Kn <- sapply(Ns, function(n) sum(colSums(data_mat[1:n,]) > 0)  )


Nbar_plot <- 400

# compute parametes of richness

# PoissonBB
list_eb_EFPF_est_PoissonBB <- list_eb_EFPF_fit_PoissonBB[grepl("Nbar.emp", names(list_eb_EFPF_fit_PoissonBB) )]
params_richness_EFPF_est_PoissonBB <- tibble( lambda_prime = unname(sapply(list_eb_EFPF_est_PoissonBB, function(x)
  total_richness(x)$lambda_post))) %>%
  add_column(n_train = Ns,
             n_train_idx = rep(c(1,2,3), each = 1),
             Model = "Poisson BB", PM = "est") %>%
  mutate(lb = qpois(0.025, lambda_prime, lower.tail = TRUE, log.p = FALSE),
         ub = qpois(0.975, lambda_prime, lower.tail = TRUE, log.p = FALSE) )

list_eb_EFPF_fixed_PoissonBB <- list_eb_EFPF_fit_PoissonBB[grepl(paste0("Nbar.",Nbar_plot), names(list_eb_EFPF_fit_PoissonBB) )]
params_richness_EFPF_fixed_PoissonBB <- tibble( lambda_prime = unname(sapply(list_eb_EFPF_fixed_PoissonBB, function(x)
  total_richness(x)$lambda_post))) %>%
  add_column(n_train = Ns,
             n_train_idx = rep(c(1,2,3), each = 1),
             Model = "Poisson BB", PM = "fixed") %>%
  mutate(lb = qpois(0.025, lambda_prime, lower.tail = TRUE, log.p = FALSE),
         ub = qpois(0.975, lambda_prime, lower.tail = TRUE, log.p = FALSE) )

# NegBinBB
params_richness_EFPF_est_NegBinBB <- tibble( n0_prime = numeric(),
                                             mu0_prime = numeric(),
                                             n_train = integer(),
                                             n_train_idx = integer(),
                                             Model = character(),
                                             PM = character(),
                                             p_prime = numeric(),
                                             lb = numeric(), ub = numeric())
params_richness_EFPF_fixed_NegBinBB <- tibble( n0_prime = numeric(),
                                               mu0_prime = numeric(),
                                               n_train = integer(),
                                               n_train_idx = integer(),
                                               Model = character(),
                                               PM = character(),
                                               p_prime = numeric(),
                                               lb = numeric(), ub = numeric())

for (var_fct_NegBinBB in vars_fct_NegBinBB){
  
  # est
  list_eb_EFPF_est_NegBinBB_var <- 
    list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]][grepl("Nbar.emp", names(list_eb_EFPF_fit_NegBinBB[[1]]) )]
  
  params_richness_EFPF_est_NegBinBB_var <- tibble( n0_prime = unname(sapply(list_eb_EFPF_est_NegBinBB_var, function(x)
    total_richness(x)$n0_post)),
    mu0_prime = unname(sapply(list_eb_EFPF_est_NegBinBB_var, function(x)
      total_richness(x)$mu0_post))) %>%
    add_column(n_train = Ns,
               n_train_idx = rep(c(1,2,3), each = 1),
               Model = paste0("NegBinomial BB x", var_fct_NegBinBB),
               PM = "est") %>%
    mutate(p_prime = 1/(mu0_prime/n0_prime + 1),
           lb = qnbinom(0.025, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE),
           ub = qnbinom(0.975, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE) )
  
  params_richness_EFPF_est_NegBinBB <- bind_rows(params_richness_EFPF_est_NegBinBB,
                                                 params_richness_EFPF_est_NegBinBB_var)
  
  
  #fixed
  list_eb_EFPF_fixed_NegBinBB_var <- 
    list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]][grepl(paste0("Nbar.",Nbar_plot), names(list_eb_EFPF_fit_NegBinBB[[1]]) )]
  
  params_richness_EFPF_fixed_NegBinBB_var <- tibble( n0_prime = unname(sapply(list_eb_EFPF_fixed_NegBinBB_var, function(x)
    total_richness(x)$n0_post)),
    mu0_prime = unname(sapply(list_eb_EFPF_fixed_NegBinBB_var, function(x)
      total_richness(x)$mu0_post))) %>%
    add_column(n_train = Ns,
               n_train_idx = rep(c(1,2,3), each = 1),
               Model = paste0("NegBinomial BB x", var_fct_NegBinBB),
               PM = "fixed") %>%
    mutate(p_prime = 1/(mu0_prime/n0_prime + 1),
           lb = qnbinom(0.025, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE),
           ub = qnbinom(0.975, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE) )
  
  
  params_richness_EFPF_fixed_NegBinBB <- bind_rows(params_richness_EFPF_fixed_NegBinBB,
                                                   params_richness_EFPF_fixed_NegBinBB_var)
  
  
  
}



bounds_distr <- tibble( lb = min(params_richness_EFPF_est_PoissonBB$lb + Kn, 
                                 params_richness_EFPF_est_NegBinBB$lb + Kn,
                                 params_richness_EFPF_fixed_PoissonBB$lb + Kn, 
                                 params_richness_EFPF_fixed_NegBinBB$lb + Kn,
                                 H),
                        ub = max(params_richness_EFPF_est_PoissonBB$ub + Kn,
                                 params_richness_EFPF_est_NegBinBB$ub + Kn,
                                 params_richness_EFPF_fixed_PoissonBB$ub + Kn,
                                 params_richness_EFPF_fixed_NegBinBB$ub + Kn,
                                 H))

# density computation
# PoissonBB
dens_richness_EFPF_est_PoissonBB <- bind_rows(lapply(1:length(Ns), function(j) 
  tibble( x = bounds_distr$lb: bounds_distr$ub) %>%
    mutate( y = dpois(x - Kn[j] , lambda = params_richness_EFPF_est_PoissonBB$lambda_prime[j])) %>%
    add_column(Model = "Poisson BB", PM = "est",
               n_train = params_richness_EFPF_est_PoissonBB$n_train[j],
               n_train_idx = params_richness_EFPF_est_PoissonBB$n_train_idx[j])
))

dens_richness_EFPF_fixed_PoissonBB <- bind_rows(lapply(1:length(Ns), function(j) 
  tibble( x = bounds_distr$lb: bounds_distr$ub) %>%
    mutate( y = dpois(x - Kn[j] , lambda = params_richness_EFPF_fixed_PoissonBB$lambda_prime[j])) %>%
    add_column(Model = "Poisson BB", PM = "fixed",
               n_train = params_richness_EFPF_fixed_PoissonBB$n_train[j],
               n_train_idx = params_richness_EFPF_fixed_PoissonBB$n_train_idx[j] )
))

# NegBinBB
dens_richness_EFPF_est_NegBinBB <- tibble(x = integer(), y = numeric(),
                                          Model = character(),
                                          PM = character(),
                                          n_train = integer(),
                                          n_train_idx = integer())

dens_richness_EFPF_fixed_NegBinBB <- tibble(x = integer(), y = numeric(),
                                            Model = character(),
                                            PM = character(),
                                            n_train = integer(),
                                            n_train_idx = integer())


for (var_fct_NegBinBB in vars_fct_NegBinBB){
  
  # est
  params_richness_EFPF_est_NegBinBB_var <- params_richness_EFPF_est_NegBinBB %>%
    filter(Model == paste0("NegBinomial BB x", var_fct_NegBinBB))
  
  dens_richness_EFPF_est_NegBinBB_var <- bind_rows(lapply(1:length(Ns), function(j) 
    tibble( x = bounds_distr$lb: bounds_distr$ub) %>%
      mutate( y = dnbinom(x - Kn[j], size = params_richness_EFPF_est_NegBinBB_var$n0_prime[j], 
                          prob = params_richness_EFPF_est_NegBinBB_var$p_prime[j])) %>%
      add_column(Model = paste0("NegBinomial BB x", var_fct_NegBinBB), 
                 PM = "est",
                 n_train = params_richness_EFPF_est_NegBinBB_var$n_train[j],
                 n_train_idx = params_richness_EFPF_est_NegBinBB_var$n_train_idx[j])
  ))
  
  dens_richness_EFPF_est_NegBinBB <- bind_rows(dens_richness_EFPF_est_NegBinBB,
                                               dens_richness_EFPF_est_NegBinBB_var)
  
  
  # fixed
  params_richness_EFPF_fixed_NegBinBB_var <- params_richness_EFPF_fixed_NegBinBB %>%
    filter(Model == paste0("NegBinomial BB x", var_fct_NegBinBB))
  
  dens_richness_EFPF_fixed_NegBinBB_var <- bind_rows(lapply(1:length(Ns), function(j) 
    tibble( x = bounds_distr$lb: bounds_distr$ub) %>%
      mutate( y = dnbinom(x - Kn[j], size = params_richness_EFPF_fixed_NegBinBB_var$n0_prime[j], 
                          prob = params_richness_EFPF_fixed_NegBinBB_var$p_prime[j])) %>%
      add_column(Model = paste0("NegBinomial BB x", var_fct_NegBinBB),
                 PM = "fixed",
                 n_train = params_richness_EFPF_fixed_NegBinBB_var$n_train[j],
                 n_train_idx = params_richness_EFPF_fixed_NegBinBB_var$n_train_idx[j])
  ))
  
  dens_richness_EFPF_fixed_NegBinBB <- bind_rows(dens_richness_EFPF_fixed_NegBinBB,
                                                 dens_richness_EFPF_fixed_NegBinBB_var)
  
}

# Set Ns and Kn to plot
Ns_plot <- Ns[1:2] # [1:2] for "custom"
Kn_plot <- Kn[1:2]

dens_richnesses <- rbind(dens_richness_EFPF_est_PoissonBB,
                         dens_richness_EFPF_est_NegBinBB,
                         dens_richness_EFPF_fixed_PoissonBB,
                         dens_richness_EFPF_fixed_NegBinBB ) %>%
  filter(n_train %in% Ns_plot)


dens_richnesses$Model <- factor(dens_richnesses$Model, 
                                levels = c("Poisson BB",
                                           paste0("NegBinomial BB x", vars_fct_NegBinBB)))

dens_richnesses$PM <- factor(dens_richnesses$PM)
levels(dens_richnesses$PM) <- c("EB", "Bayesian")
dens_richnesses <- dens_richnesses %>%
  mutate(n_train_latex = paste("Scenario~", LETTERS[n_train_idx], "~~(n == ", n_train, "~~ K[n] == ", Kn[n_train_idx], ")", sep = ""))

ggplot(dens_richnesses, aes(x = x, y = y, color = Model)) +
  geom_line() +
  geom_vline(aes(xintercept = H), linetype="dashed") +
  facet_grid(PM ~ n_train_latex,
             labeller = label_parsed,
             scales = "free") +
  theme_light() +
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  xlab("# distinct features") + ylab("Probability") + 
  scale_color_tableau()
