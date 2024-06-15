rm(list=ls())

library(ProductFormFA)
library(tidyverse)
library(scales)
library(ggpubr)
library(ggthemes)

source("R_script_paper/Routine_Chao.R")
source("R_script_paper/utils.R")


# Choose dataset among the 3 available
type = "Plants" # options: "Lichens", "Plants" 

data <- read.csv(file = paste0("R_script_paper/mazz2016_",type,".csv"), header = TRUE,
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
  geom_point(color="black", shape = 21, size = 1) + 
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) + 
  scale_color_tableau()
ggsave(filename = paste0("R_script_paper/Paper_plots/accumulation_", type, ".pdf"), width = 4, height = 4, dpi = 300, units = "in", device='pdf')


vars_fct_NegBinBB <- c(10, 1000)
vars_GammaIBP <- c(0.01, 100)

# EFPF approach: fit models -----

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


# Model Checking ----

# 0.B) Check on rarefaction
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

# for plot
rare_EFPF_mixtureBB <- rare_EFPF_NegBinBB %>%
  mutate(Model = "PoissonBB/NegBinBB")

df_rare <- rbind(rare_EFPF_mixtureBB,
                 rare_EFPF_GammaIBP)

df_rare$Model <- factor(df_rare$Model,
                        levels = c("PoissonBB/NegBinBB", "GammaIBP"))

ggplot(accum_df, aes(x = x, y = n_feat)) +
  geom_point(color="black", shape = 18, size = 0.8) + 
  #facet_wrap(.~ Model, scales = "free_x", nrow = 1) +
  #geom_ribbon(aes(ymin = lb_bands, ymax = ub_bands), color = "red" , linewidth = 0.8, alpha = 0.1) +
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
ggsave(filename = paste0("R_script_paper/Paper_plots/rarefaction_", type, "_eb_EFPF.pdf"), width = 4, height = 4, dpi = 300, units = "in", device='pdf')


# 0.c) Check on K_n_r
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

# for plot
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
ggsave(filename = paste0("R_script_paper/Paper_plots/knr_", type, "_eb_EFPF.pdf"), width = 4, height = 4, dpi = 300, units = "in", device='pdf')



# # Prediction ----
# 
# ## Richness estimation -----
# 
# params_richness_EFPF_PoissonBB <- tibble( lambda_prime = 
#                                             total_richness(eb_EFPF_PoissonBB_fit)$lambda_post) %>%
#   add_column(Model = "PoissonBB") %>%
#   mutate(lb = qpois(0.025, lambda_prime, lower.tail = TRUE, log.p = FALSE),
#          ub = qpois(0.975, lambda_prime, lower.tail = TRUE, log.p = FALSE) )
# 
# 
# params_richness_EFPF_NegBinBB <- tibble( n0_prime = 
#                                            total_richness(eb_EFPF_NegBinBB_fit)$n0_post,
#                                          mu0_prime = 
#                                            total_richness(eb_EFPF_NegBinBB_fit)$mu0_post) %>%
#   add_column(Model = "NegBinBB") %>%
#   mutate(p_prime = 1/(mu0_prime/n0_prime + 1),
#          lb = qnbinom(0.025, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE),
#          ub = qnbinom(0.975, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE) )
# 
# 
# 
# bounds <- tibble( lb = min(params_richness_EFPF_PoissonBB$lb + Kn, 
#                            params_richness_EFPF_NegBinBB$lb + Kn),
#                   ub = max(params_richness_EFPF_PoissonBB$ub + Kn,
#                            params_richness_EFPF_NegBinBB$ub + Kn))
# 
# 
# dens_richness_PoissonBB <- tibble( x = bounds$lb: bounds$ub) %>%
#   mutate( y = dpois(x - Kn , lambda = params_richness_EFPF_PoissonBB$lambda_prime)) %>%
#   add_column(Model = "PoissonBB")
# 
# dens_richness_NegBinBB <- tibble( x = bounds$lb : bounds$ub) %>%
#   mutate( y = dnbinom(x - Kn, size = params_richness_EFPF_NegBinBB$n0_prime, 
#                       prob = params_richness_EFPF_NegBinBB$p_prime)) %>%
#   add_column(Model = "NegBinBB")
# 
# dens_richnesses <- rbind(dens_richness_PoissonBB,
#                          dens_richness_NegBinBB)
# 
# 
# ggplot(dens_richnesses, aes(x = x, y = y, color = Model)) +
#   geom_line() +
#   theme_light() +
#   facet_wrap(~"Richness") +
#   theme(legend.position = "top") +
#   scale_y_continuous(breaks = pretty_breaks()) +
#   xlab("# distinct features") + rremove("ylab") +
#   scale_color_tableau() +
#   theme(aspect.ratio = 1)
# ggsave(filename = paste0("R_script_paper/Paper_plots/richness_", type, "_eb_EFPF.pdf"), width = 3.8, height = 3.8, dpi = 300, units = "in", device='pdf')


## Extrapolation  -----

# alpha-diversity of the gamma mixture
idx_gamma <- 2

a <- list_eb_EFPF_fit_GammaIBP[[idx_gamma]]$a
b <- list_eb_EFPF_fit_GammaIBP[[idx_gamma]]$b
alpha <- list_eb_EFPF_fit_GammaIBP[[idx_gamma]]$alpha
theta <- list_eb_EFPF_fit_GammaIBP[[idx_gamma]]$theta

gamma_a_t_n <- sum(exp(lgamma(alpha + theta + (1:n) - 1) - lgamma(alpha + theta) -
                         lgamma(theta + (1:n)) + lgamma(theta +1) ) )

a_gamma <- a + Kn
b_gamma <- (b + gamma_a_t_n)*gamma(theta+alpha)/gamma(theta+1)*alpha

print(paste0("a_gamma = ", a_gamma))  
print(paste0("b_gamma = ", b_gamma))  

print(paste0("mean = ", a_gamma/b_gamma))  
print(paste0("var = ", a_gamma/b_gamma^2))  

# Extract accumulation curve of the observed sample (or average accumulation)
M = 400

accum_df <- tibble( x = 0:n,
                    n_feat = c(0,rarefaction(data_mat, n_reorderings = 200)))


# # PoissonBB
# extr_EFPF_PoissonBB_df <- tibble(lambda = unname(unlist( 
#   extrapolation(object = eb_EFPF_fit_PoissonBB, M = M, seed = seed)$lambda_post)),
#   Kn = rep(Kn, each = M)) %>%
#   mutate(lb = qpois(0.025, lambda, lower.tail = TRUE, log.p = FALSE),
#          ub = qpois(0.975, lambda, lower.tail = TRUE, log.p = FALSE)) %>%
#   rename(means = lambda) %>%
#   add_row(means = 0, lb = 0, ub = 0, Kn = Kn) %>%
#   mutate(means = means + Kn, lb = lb + Kn, ub = ub + Kn) %>%
#   add_column(x = c((n+1):(n+M), n),
#              Model = "PoissonBB")
# 
# extr_EFPF_PoissonBB_df$x <- as.integer(extr_EFPF_PoissonBB_df$x)
# extr_EFPF_PoissonBB_df <- extr_EFPF_PoissonBB_df %>%
#   select(means, lb, ub, x, Model)
# 
# # NegBin
# extr_EFPF_NegBinBB_df <- tibble(means = numeric(), 
#                                 lb = numeric(), ub = numeric(),
#                                 x = integer(), Model = character())
# 
# for (var_fct_NegBinBB in vars_fct_NegBinBB){
#   
#   eb_EFPF_NegBinBB_var <- list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]]
#   
#   extr_EFPF_NegBinBB_df_var <- tibble(mu0 = unname(unlist(
#     extrapolation(object = eb_EFPF_NegBinBB_var, M = M, seed = seed)$mu0_post )),
#     n0 = unname(unlist( extrapolation(object = eb_EFPF_NegBinBB_var, M = M, seed = seed)$n0_post )),
#     Kn = rep(Kn, each = M)) %>%
#     mutate(p = 1/(mu0/n0 + 1),
#            lb = qnbinom(0.025, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE),
#            ub = qnbinom(0.975, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE)) %>%
#     rename(means = mu0) %>%
#     add_row(means = 0, lb = 0, ub = 0, Kn = Kn) %>%
#     mutate(means = means + Kn, lb = lb + Kn, ub = ub + Kn) %>%
#     add_column(x = c((n+1):(M+n), n),
#                Model = paste0("NegBinBB x", var_fct_NegBinBB))
#   
#   extr_EFPF_NegBinBB_df_var$x <- as.integer(extr_EFPF_NegBinBB_df_var$x)
#   extr_EFPF_NegBinBB_df_var <- extr_EFPF_NegBinBB_df_var %>%
#     select(means, lb, ub, x, Model)
#   
#   extr_EFPF_NegBinBB_df <- bind_rows(extr_EFPF_NegBinBB_df, 
#                                      extr_EFPF_NegBinBB_df_var)
# }

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

# extr_EFPF_PoissonBB_df <- extr_EFPF_PoissonBB_df %>%
#   add_column(Model_gen = "PoissonBB/NegBinBB")
# extr_EFPF_NegBinBB_df <- extr_EFPF_NegBinBB_df %>%
#   add_column(Model_gen = "PoissonBB/NegBinBB")
extr_EFPF_GammaIBP_df <- extr_EFPF_GammaIBP_df %>%
  add_column(Model_gen = "GammaIBP")
extr_all_df <- rbind(#extr_EFPF_PoissonBB_df, 
                      #extr_EFPF_NegBinBB_df,
                     extr_EFPF_GammaIBP_df)


extr_all_df$Model <- factor(extr_all_df$Model,
                            levels = c("Poisson BB",
                                       paste0("NegBinomial BB x", vars_fct_NegBinBB),
                                       paste0("Gamma IBP, Variance: ", vars_GammaIBP)))

extr_all_df$Model_gen <- factor(extr_all_df$Model_gen,
                                levels = c("PoissonBB/NegBinBB", "GammaIBP"))

ggplot(extr_all_df, aes(x, means, color = Model)) +
  geom_line(linetype = "dashed") +
  #facet_wrap(. ~ Model_gen,  scales = "free_x") +
  geom_point( data = accum_df, aes(x = x, y = n_feat),
              color="black", shape = 19, size = 0.1) +
  geom_ribbon(aes(ymin = lb, ymax = ub, color = Model), alpha = 0) +
  #geom_line(data = df_extr_GT_long, aes(t, value)) +
  #geom_line(data = df_extr_Chao_long, aes(t, medians), linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = n) , linetype = "dashed", color = "grey") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  #facet_wrap(~"Extrapolation") +
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) +
  scale_color_tableau()
ggsave(filename = paste0("R_script_paper/Paper_plots/extr_", type, "_eb_EFPF.pdf"), width = 5, height = 5, dpi = 300, units = "in", device='pdf')




extr_EFPF_GammaIBP_df %>%
  filter(Model == "GammaIBP, var = 0.01",
         x %in% c(n + 1, n + 10, n + 100, n + 1000)) %>%
  mutate(means_new = means - Kn,
         lb_new = lb -Kn,
         ub_new = ub - Kn)



# Prior approach ----------

# We focus on GammaIBP + prior (since it is selected from model-checking)
  
# Fit for GibbsFA's (save workspace)
if (!file.exists(paste0("R_script_paper/prior_",type,"_fit_estimate_singledataset.RData"))) {
  
  list_prior_fit_GammaIBP <-  vector(mode = "list", length = length(vars_GammaIBP))
  names(list_prior_fit_GammaIBP) <- paste0("var.", vars_GammaIBP)
  
  # Initialization and MCMC setting 
  mcmcparams_GammaIBP <- list(sigq_alpha = 0.1, sigq_s = 0.1, 
                              S = 5*10^4, n_burnin = 5*10^3, thin = 2)
  mcmcparams_obj_GammaIBP <- mcmcparameters(model = "GammaIBP", mcmcparams = mcmcparams_GammaIBP)
  
  # 1) more_prior
  # init_GammaIBP <- list(alpha_0 = 0.5, s_0 = 15, a_0 = 5, b_0 = 1)
  # init_obj_GammaIBP <- initialization(model = "GammaIBP_more_prior", init = init_GammaIBP )
  
  # 2) single_prior
  init_GammaIBP <- list(alpha_0 = 0.5, s_0 = 15)
  init_obj_GammaIBP <- initialization(model = "GammaIBP_single_prior", init = init_GammaIBP )
  
  # EB estimates
  small_val <- 2
  alpha_eb <- list_eb_EFPF_fit_GammaIBP[[1]]$alpha
  theta_eb <- list_eb_EFPF_fit_GammaIBP[[1]]$theta
  
  t_eb <- (1 - alpha_eb)/alpha_eb
  s_eb <- alpha_eb + theta_eb
  
  print(paste0("Prior variance of alpha: ", 
               t_eb/(1 + t_eb)^2 /(1 + small_val*(1+t_eb))))
  
  # Fit the model
  for (var_GammaIBP in vars_GammaIBP){
    
    a_eb <- list_eb_EFPF_fit_GammaIBP[[paste0("var.", var_GammaIBP)]]$a
    b_eb <- list_eb_EFPF_fit_GammaIBP[[paste0("var.", var_GammaIBP)]]$b
    
    # Hyperparameters elicitation 
    # 1) more_prior
    # hyper_GammaIBP <- list(a_alpha = small_val, b_alpha = t_eb*small_val,
    #                        a_s = s_eb*small_val , b_s = small_val,
    #                        q = 1/a_eb, r = b_eb*small_val, t = small_val)
    # prior_obj_GammaIBP <- prior(model = "GammaIBP_more_prior", hyper = hyper_GammaIBP)
    # 
    # 2) single_prior
    hyper_GammaIBP <- list(a = a_eb, b = b_eb,
                           a_alpha = small_val, b_alpha = t_eb*small_val,
                           a_s = s_eb*small_val , b_s = small_val)
    prior_obj_GammaIBP <- prior(model = "GammaIBP_single_prior", hyper = hyper_GammaIBP)
    
    
    list_prior_fit_GammaIBP[[paste0("var.", var_GammaIBP)]] <- 
      GibbsFA(feature_matrix = data_mat,
              model = "GammaIBP_single_prior", 
              prior = prior_obj_GammaIBP,
              initialization = init_obj_GammaIBP,
              mcmcparams = mcmcparams_obj_GammaIBP)
    
  }
  
  # Save the entire workspace related to the type just performed
  save(list = ls(all.names = TRUE), file =  paste0("R_script_paper/prior_",type,"_fit_estimate_singledataset.RData"))
  
}


# Load the Work space
load(paste0("R_script_paper/prior_",type,"_fit_estimate_singledataset.RData"))


##### Convergence checks -------------
library(ggmcmc)
library(coda)

params_prior_GammaIBP <- list_prior_fit_GammaIBP[[paste0("var.", vars_GammaIBP[2])]][c("a_chain","b_chain", "alpha_chain", "theta_chain")]
params_prior_GammaIBP_df <- as.data.frame(do.call(cbind, params_prior_GammaIBP))

samples_GammaIBP <- mcmc.list(mcmc(params_prior_GammaIBP_df))
samples_ggs_GammaIBP <- ggs(samples_GammaIBP, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_GammaIBP)

effectiveSize(params_prior_GammaIBP_df)


##### Prior: Extrapolation ------

if (!file.exists(paste0("R_script_paper/prior_",type,"_extrapolation.RData"))) {
  
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
  save(list = ls(all.names = TRUE), file =  paste0("R_script_paper/prior_",type,"_extrapolation.RData"))
  
}

# Load the Work space
load(paste0("R_script_paper/prior_",type,"_extrapolation.RData"))
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


ggplot(extr_joint_GammaIBP_df, aes(x, means, color = Type )) +
  geom_line(linetype = "dashed") +
  facet_wrap(. ~ Model,  scales = "free_x") +
  geom_point( data = accum_df, aes(x = x, y = n_feat),
              color="black", shape = 19, size = 0.1) +
  geom_ribbon(aes(ymin = lb, ymax = ub, color = Type), alpha = 0) +
  #geom_line(data = df_extr_GT_long, aes(t, value)) +
  #geom_line(data = df_extr_Chao_long, aes(t, medians), linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = n) , linetype = "dashed", color = "grey") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  labs(color = "Approach") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) +
  scale_color_tableau()
ggsave(filename = paste0("R_script_paper/Paper_plots/extr_", type, "_prior.pdf"), width = 6, height = 4, dpi = 300, units = "in", device='pdf')



### OLD STUFF -----

# 3) Estimate the total richness ----

PoissonBB_rich <- total_richness(object = PoissonBB_fit)
eb_PoissonBB_rich <- total_richness(object = eb_PoissonBB_fit)
plot(x = eb_PoissonBB_fit, type = "richness")

NegBinBB_rich <- total_richness(object = NegBinBB_fit)
eb_NegBinBB_rich <- total_richness(object = eb_NegBinBB_fit)
plot(x = eb_NegBinBB_fit, type = "richness")

# 4) Extrapolation ----

PoissonBB_extr <- extrapolation(object = PoissonBB_fit, M = M) 
NegBinBB_extr <- extrapolation(object = NegBinBB_fit, M = M) 

# 5) Rarefaction and checks ----
emp_pis <- data.frame(x = colMeans(data_mat))

ggplot(emp_pis, aes(x = x) ) +
  stat_ecdf(linewidth=2, colour = "red") +
  stat_function(fun = pbeta, colour = "blue", 
                args = list(shape1 = - eb_PoissonBB_fit$alpha,
                            shape2 = eb_PoissonBB_fit$alpha + eb_PoissonBB_fit$theta)) +
  labs(title="ECDF and theoretical CDF")  


PoissonBB_rare <- rarefaction(object = PoissonBB_fit) 
eb_PoissonBB_rare <- rarefaction(object = eb_PoissonBB_fit) 
plot(x = eb_PoissonBB_fit, type = "rarefaction", n_reorderings = 10)

NegBinBB_rare <- rarefaction(object = NegBinBB_fit) 
eb_NegBinBB_rare <- rarefaction(object = eb_NegBinBB_fit)
plot(x = eb_NegBinBB_fit, type = "rarefaction", n_reorderings = 10)

# 6) Competitors: Chao and GT ------

# 6.1) Chao: competitor for rarefaction and extrapolation

# Determine the frequency vector of the training sets
Q_vec <- colSums(data_mat)
Q_vec <- Q_vec[Q_vec>0]

# Compute the curves with confidence intervals
fit_Chao <- iNEXT.Sam(Spec = Q_vec, T = n, endpoint = n + M)

rare_Chao <- as_tibble(fit_Chao[["q=0"]]) %>%
  select(-Cov.hat) %>%
  rename(medians = D0.hat, lbs = Norm.CI.Low, ubs = Norm.CI.High)

Chao_rare_extr <- as.data.frame(rare_Chao)


# 6.2) Smoothed Good-Toulmin: competitor for extrapolation

# Compute SFS vector and CTS vector
sfs <- tabulate(colSums(data_mat))
cts <- sapply(2:n, function(i) ncol(data_mat[1:i,colSums(data_mat[1:i,]) > 0])   )
cts <- c(0, sum(data_mat[1,]) , cts)

GT_extr <- predict_good_toulmin(n, M, sfs, cts, alternative = 0)$preds



# 7) Plot Richness: whole distributions for Poisson and NegBin

richness_PoissonBB_long <- as_tibble(PoissonBB_rich) %>%
  add_column(model= "Poisson")

richness_NegBinBB_long <- as_tibble(NegBinBB_rich) %>%
  add_column(model= "NegBin")

joint_richness_long <- bind_rows(richness_PoissonBB_long, richness_NegBinBB_long) %>%
  mutate(model = fct_relevel(model, c( "NegBin", "Poisson")))


ggplot(joint_richness_long, aes(x = value, color = model)) +
  stat_density(aes(x=value, colour=model),
               geom="line",position="identity", adjust = 2.5) +
  theme_light() +
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  xlab("# distinct features") + rremove("ylab") +
  scale_color_tableau() +
  theme(aspect.ratio = 1)



# 8) Plot Rarefaction: PoissonBB, NegBinBB and Chao

# PoissonBB
df_rare_PoissonBB <- as_tibble(t(bind_rows(as.data.frame(lapply(PoissonBB_rare, quantile, prob = c(0.025, 0.975))),
                                           as.data.frame(lapply(PoissonBB_rare, mean))))) 
colnames(df_rare_PoissonBB) <- c("lbs", "ubs", "means")
df_rare_PoissonBB <- df_rare_PoissonBB %>%
  add_column(t = 1:nrow(df_rare_PoissonBB),
             model = "Poisson") %>%
  add_row(means = 0, ubs = 0, lbs = 0, t=0, model = "Poisson")


# NegBinBB
df_rare_NegBinBB <- as_tibble(t(bind_rows(as.data.frame(lapply(NegBinBB_rare, quantile, prob = c(0.025, 0.975))),
                                          as.data.frame(lapply(NegBinBB_rare, mean))))) 
colnames(df_rare_NegBinBB) <- c("lbs", "ubs", "means")
df_rare_NegBinBB <- df_rare_NegBinBB %>%
  add_column(t = 1:nrow(df_rare_NegBinBB),
             model = "NegBin") %>%
  add_row(means = 0, ubs = 0, lbs = 0, t=0, model = "NegBin")


# Chao
df_rare_Chao <- Chao_rare_extr %>%
  select(-c(lbs,ubs)) %>%
  rename(value = medians) %>%
  add_column(model = "Chao") %>%
  add_row(value = 0, t=0, model = "Chao") %>%
  filter(t <= n)


# Accumulation curve
accum <- rarefaction(data_mat)
accum_df <- data.frame("accum" = c(0, accum),
                       "t" = 0:length(accum))

joint_df_rare_bayes <- rbind(df_rare_PoissonBB,df_rare_NegBinBB)

# plot
ggplot(joint_df_rare_bayes, aes(x = t, y = means, color = model)) +
  geom_line(linetype = "dashed", linewidth = 0.8) +
  geom_ribbon(aes(ymin = lbs, ymax = ubs), linewidth = 0.8, alpha = 0.1) +
  geom_line(data = df_rare_Chao, aes(t, value), linewidth = 0.8) +
  geom_line(data = accum_df, aes(t, accum), color="black", linetype="solid") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_tableau() +
  theme(aspect.ratio = 1)




# 9) Plot Extrapolation: PoissonBB, NegBinBB, Chao and GT

# Set extrapolation horizon 
M <- 1000

# PoissonBB
df_extr_PoissonBB <- as_tibble(t(bind_rows(as.data.frame(lapply(PoissonBB_extr, quantile, prob = c(0.025, 0.975))),
                                           as.data.frame(lapply(PoissonBB_extr, mean))))) 
colnames(df_extr_PoissonBB) <- c("lbs", "ubs", "means")
df_extr_PoissonBB <- df_extr_PoissonBB %>%
  add_column(t = 1:nrow(df_extr_PoissonBB),
             model = "Poisson") %>%
  mutate( t = t + n ) %>%
  add_row(means = Kn, ubs = Kn, lbs = Kn, t=n, model = "Poisson")


# NegBinBB
df_extr_NegBinBB <- as_tibble(t(bind_rows(as.data.frame(lapply(NegBinBB_extr, quantile, prob = c(0.025, 0.975))),
                                          as.data.frame(lapply(NegBinBB_extr, mean))))) 
colnames(df_extr_NegBinBB) <- c("lbs", "ubs", "means")
df_extr_NegBinBB <- df_extr_NegBinBB %>%
  add_column(t = 1:nrow(df_extr_NegBinBB),
             model = "NegBin") %>%
  mutate( t = t + n ) %>%
  add_row(means = Kn, ubs = Kn, lbs = Kn, t=n, model = "NegBin")


# Chao
df_extr_Chao <- Chao_rare_extr %>%
  select(-c(lbs,ubs)) %>%
  rename(value = medians) %>%
  add_column(model = "Chao") %>%
  filter(t >= n)

# Good-Toulmin
df_extr_GT <- as_tibble(GT_extr) %>%
  add_column(t = 0:(length(GT_extr)-1),
             model = "GT") %>%
  filter(t >= n)


# Accumulation curve
accum <- rarefaction(data_mat)
accum_df <- data.frame("accum" = c(0, accum),
                       "t" = 0:length(accum))


joint_df_extr_bayes <- rbind(df_extr_PoissonBB,df_extr_NegBinBB)

# plot
ggplot(joint_df_extr_bayes, aes(x = t, y = means, color = model)) +
  geom_line(linetype = "dashed", linewidth = 0.8) +
  geom_ribbon(aes(ymin = lbs, ymax = ubs), linewidth = 0.8, alpha = 0.1) +
  geom_line(data = accum_df, aes(t, accum), color="black", linetype="solid") +
  geom_line(data = df_extr_Chao, aes(t, value), linewidth = 0.8) +
  geom_line(data = df_extr_GT, aes(t, value), linewidth = 0.8) +
  geom_vline(aes(xintercept = n), linetype="dashed", color = "grey") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_tableau() +
  theme(aspect.ratio = 1)








