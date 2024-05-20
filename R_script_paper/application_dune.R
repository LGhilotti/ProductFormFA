rm(list=ls())

library(ProductFormFA)
library(tidyverse)
library(scales)
library(ggpubr)
library(ggthemes)

source("R_script_paper/Routine_Chao.R")
source("R_script_paper/utils.R")


library(vegan)
data(dune)

data <- 1*as.matrix(dune > 0)

data <- data[, colSums(is.na(data))==0]
data <- data[, colSums(data)!=0]


# Number of sites and number of species
n <- nrow(data)
Kn <- ncol(data)
print(paste0("Number of sites: ", n ))
print(paste0("Number of species: ", Kn))

# Set extrapolation horizon 
M <- 1000

# Randomly reorder sites
seed <- 12345
set.seed(seed)

data_mat <- data[sample.int(n, size = n, replace = F),]

# Plot accumulation
plot(1:n, rarefaction(data_mat, n_reorderings = 10) )


# 2) Run the models on the data ----

# Empirical estimate of E(N) is obtained by Chiu
Nbar <- beta_binomial_estimator(data_mat)


# 2.1) Beta-Bernoulli with Poisson(lambda) mixture - EFPF version

# EB parameters 
eb_init_PoissonBB <- list(alpha = -100, s = 10, lambda = 100)
eb_known_PoissonBB <- list()
eb_params_obj_PoissonBB <- eb_params(model = "PoissonBB", 
                                     init = eb_init_PoissonBB, known = eb_known_PoissonBB )

eb_PoissonBB_fit <- GibbsFA_eb(feature_matrix = data_mat, 
                               model = "PoissonBB", 
                               type = "EFPF",
                               eb_params =  eb_params_obj_PoissonBB)


eb_PoissonBB_fit[c("alpha", "theta", "lambda")]

# 2.1:B) Beta-Bernoulli with Poisson(lambda) mixture - MM_biased version

# EB parameters 
eb_PoissonBB_MM_biased_fit <- GibbsFA_eb(feature_matrix = data_mat, 
                               model = "PoissonBB", 
                               type = "MM_biased")


eb_PoissonBB_MM_biased_fit[c("alpha", "theta", "lambda")]
- eb_PoissonBB_MM_biased_fit[["alpha"]]/eb_PoissonBB_MM_biased_fit[["theta"]]


# 2.1:B) Beta-Bernoulli with Poisson(lambda) mixture - MM_censored version

# EB parameters 
eb_PoissonBB_MM_cens_fit <- GibbsFA_eb(feature_matrix = data_mat, 
                                  model = "PoissonBB", 
                                  type = "MM_censored")


eb_PoissonBB_MM_cens_fit[c("alpha", "theta", "lambda")]
- eb_PoissonBB_MM_cens_fit[["alpha"]]/eb_PoissonBB_MM_cens_fit[["theta"]]

# 2.2) NegBinBB - EFPF version

# Initialization and known parameters
c_fr <- 10
eb_init_NegBinBB <- list(alpha = -1, s = 4, mu0 = 3)
eb_known_NegBinBB <- list(var_fct = c_fr)
eb_params_obj_NegBinBB <- eb_params(model = "NegBinBB",
                                    init = eb_init_NegBinBB,
                                    known = eb_known_NegBinBB)

eb_NegBinBB_fit <- GibbsFA_eb(feature_matrix = data_mat,
                              model = "NegBinBB", type = "EFPF",
                              eb_params =  eb_params_obj_NegBinBB)

eb_NegBinBB_fit[c("alpha", "theta", "mu0", "n0")]

# 2.2.B) NegBinBB - MM version
var_fct <- c_fr
eb_NegBinBB_MM_fit <- GibbsFA_eb(feature_matrix = data_mat,
                              model = "NegBinBB",
                              type = "MM", var_fct)

eb_NegBinBB_MM_fit[c("alpha", "theta", "mu0", "n0")]


# Different variances
cfrs <- c(2, 10, 200, 1000)
eb_NegBinBB_list_variances <- vector(mode="list", length = length(cfrs))
names(eb_NegBinBB_list_variances) <- paste0("c_fr.", cfrs)

eb_NegBinBB_MM_list_variances <- vector(mode="list", length = length(cfrs))
names(eb_NegBinBB_MM_list_variances) <- paste0("c_fr.", cfrs)

for (cfr in cfrs){
  eb_known_NegBinBB_vars <- list(var_fct = cfr)
  eb_params_obj_NegBinBB <- eb_params(model = "NegBinBB",
                                      init = eb_init_NegBinBB,
                                      known = eb_known_NegBinBB_vars)
  
  eb_NegBinBB_list_variances[[paste0("c_fr.", cfr)]] <- GibbsFA_eb(feature_matrix = data_mat,
                                                         model = "NegBinBB", type = "EFPF",
                                                         eb_params =  eb_params_obj_NegBinBB)
  
  var_fct <- cfr
  eb_NegBinBB_MM_list_variances[[paste0("c_fr.", cfr)]] <- GibbsFA_eb(feature_matrix = data_mat,
                                                                      model = "NegBinBB", type = "MM",
                                                                      var_fct)
}


# 2.3) GammaIBP

# Initialization and known parameters
eb_init_GammaIBP <- list(alpha = 0.5, s = 1, a = 1, b = 1)
eb_known_GammaIBP <- list()
eb_params_obj_GammaIBP <- eb_params(model = "GammaIBP",
                                    init = eb_init_GammaIBP,
                                    known = eb_known_GammaIBP)

# EFPF based
eb_GammaIBP_EFPF_fit <- GibbsFA_eb(feature_matrix = data_mat,
                              model = "GammaIBP", type = "EFPF",
                              eb_params =  eb_params_obj_GammaIBP)

eb_GammaIBP_EFPF_fit[c("alpha", "theta", "a", "b")]
eb_GammaIBP_EFPF_fit[["a"]]/eb_GammaIBP_EFPF_fit[["b"]]
eb_GammaIBP_EFPF_fit[["a"]]/eb_GammaIBP_EFPF_fit[["b"]]^2


# MM based
eb_GammaIBP_MM_fit <- GibbsFA_eb(feature_matrix = data_mat,
                                   model = "GammaIBP", type = "MM",
                                   eb_params =  eb_params_obj_GammaIBP)

eb_GammaIBP_MM_fit[c("alpha", "theta", "a", "b")]
eb_GammaIBP_MM_fit[["a"]]/eb_GammaIBP_MM_fit[["b"]]
eb_GammaIBP_MM_fit[["a"]]/eb_GammaIBP_MM_fit[["b"]]^2

# just check that the function evaluation of mse_GammaIBP is higher in EFPF estimates
eb_init_GammaIBP_eval <- list()
eb_known_GammaIBP_eval <- list(alpha = eb_GammaIBP_EFPF_fit[["alpha"]], 
                          s = eb_GammaIBP_EFPF_fit[["alpha"]] + eb_GammaIBP_EFPF_fit[["theta"]],
                          a = eb_GammaIBP_EFPF_fit[["a"]], 
                          b = eb_GammaIBP_EFPF_fit[["b"]])

eb_params_obj_GammaIBP_eval <- eb_params(model = "GammaIBP",
                                    init = eb_init_GammaIBP_eval,
                                    known = eb_known_GammaIBP_eval)

eb_GammaIBP_MM_evaluation <- GibbsFA_eb(feature_matrix = data_mat,
                                 model = "GammaIBP", type = "MM",
                                 eb_params =  eb_params_obj_GammaIBP_eval)

eb_GammaIBP_MM_evaluation$fun_value
eb_GammaIBP_MM_fit$fun_value

# Model Checking

emp_pis <- data.frame(x = colMeans(data_mat))

# PoissonBB
grid <- seq(0,1, length.out = 1000)

a_beta_PoissonBB <- - eb_PoissonBB_fit$alpha
b_beta_PoissonBB <- eb_PoissonBB_fit$alpha + eb_PoissonBB_fit$theta

cdf_betas_cond_PoissonBB <- data.frame(
  x = grid,
  y = pbeta(grid, shape1 = a_beta_PoissonBB, shape2 = b_beta_PoissonBB)*
    (1/(1- beta(a_beta_PoissonBB, n + b_beta_PoissonBB)/beta(a_beta_PoissonBB, b_beta_PoissonBB))) +
    pbeta(grid, shape1 = a_beta_PoissonBB, shape2 = n + b_beta_PoissonBB)*
    (1/(1- beta(a_beta_PoissonBB, b_beta_PoissonBB)/beta(a_beta_PoissonBB, n + b_beta_PoissonBB)))
  )

# cdf_beta <- data.frame(x = grid,
#                        y = pbeta(grid, shape1 = - eb_PoissonBB_fit$alpha,
#                                  shape2 = eb_PoissonBB_fit$alpha + eb_PoissonBB_fit$theta))

ggplot(emp_pis, aes(x = x) ) +
  stat_ecdf(linewidth=2, colour = "red") +
  geom_line(data = cdf_betas_cond_PoissonBB, aes(x = x, y = y), colour = "blue") +
  labs(title="ECDF and theoretical CDF")  


plot(x = eb_PoissonBB_fit, type = "richness")


# NegBinBB
grid <- seq(0,1, length.out = 1000)

a_beta_NegBinBB <- - eb_NegBinBB_fit$alpha
b_beta_NegBinBB <- eb_NegBinBB_fit$alpha + eb_NegBinBB_fit$theta

cdf_betas_cond_NegBinBB <- data.frame(
  x = grid,
  y = pbeta(grid, shape1 = a_beta_NegBinBB, shape2 = b_beta_NegBinBB)*
    (1/(1- beta(a_beta_NegBinBB, n + b_beta_NegBinBB)/beta(a_beta_NegBinBB, b_beta_NegBinBB))) +
    pbeta(grid, shape1 = a_beta_NegBinBB, shape2 = n + b_beta_NegBinBB)*
    (1/(1- beta(a_beta_NegBinBB, b_beta_NegBinBB)/beta(a_beta_NegBinBB, n + b_beta_NegBinBB)))
)
# 
# cdf_beta <- data.frame(x = grid,
#                        y = pbeta(grid, shape1 = - eb_NegBinBB_fit$alpha,
#                                  shape2 = eb_NegBinBB_fit$alpha + eb_NegBinBB_fit$theta))

ggplot(emp_pis, aes(x = x) ) +
  stat_ecdf(linewidth=2, colour = "red") +
  geom_line(data = cdf_betas_cond_NegBinBB, aes(x = x, y = y), colour = "blue") +
  labs(title="ECDF and theoretical CDF")  


plot(x = eb_NegBinBB_fit, type = "richness")


###### Check on rarefaction -------------

# 3 models 

accum_df <- tibble( x = 0:n,
                    n_feat = c(0,rarefaction(data_mat[1:n,], n_reorderings = 20)))

rare_PoissonBB_EFPF <- tibble( lambda_post = unname(unlist(
  rarefaction(object = eb_PoissonBB_fit, seed = seed)$lambda_post ))) %>%
  rename(means = lambda_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "PoissonBB_EFPF",
             x = c(1:n,0))

rare_PoissonBB_MM_biased <- tibble( lambda_post = unname(unlist(
  rarefaction(object = eb_PoissonBB_MM_biased_fit, seed = seed)$lambda_post ))) %>%
  rename(means = lambda_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "PoissonBB_MM_biased",
             x = c(1:n,0))

rare_PoissonBB_MM_cens <- tibble( lambda_post = unname(unlist(
  rarefaction(object = eb_PoissonBB_MM_cens_fit, seed = seed)$lambda_post ))) %>%
  rename(means = lambda_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "PoissonBB_MM_cens",
             x = c(1:n,0))

# rare_NegBinBB_EFPF <- tibble( mu0_post = unname(unlist(
#   rarefaction(object = eb_NegBinBB_fit, seed = seed)$mu0_post ))) %>%
#   rename(means = mu0_post) %>%
#   add_row(means = 0) %>%
#   add_column(Model = "NegBinBB_EFPF",
#              x = c(1:n,0))
# 
# rare_NegBinBB_MM <- tibble( mu0_post = unname(unlist(
#   rarefaction(object = eb_NegBinBB_MM_fit, seed = seed)$mu0_post ))) %>%
#   rename(means = mu0_post) %>%
#   add_row(means = 0) %>%
#   add_column(Model = "NegBinBB_MM",
#              x = c(1:n,0))
# 

rare_GammaIBP_EFPF <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_GammaIBP_EFPF_fit, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "GammaIBP_EFPF",
             x = c(1:n,0))

rare_GammaIBP_MM <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_GammaIBP_MM_fit, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "GammaIBP_MM",
             x = c(1:n,0))

df_rare <- rbind(rare_PoissonBB_EFPF, rare_PoissonBB_MM_biased, rare_PoissonBB_MM_cens,
                 #rare_NegBinBB_EFPF, rare_NegBinBB_MM, 
                 rare_GammaIBP_EFPF, rare_GammaIBP_MM)

df_rare$Model <- factor(df_rare$Model)

ggplot(df_rare, aes(x = x, y = means, color = Model)) +
  geom_line(linetype = "solid", color = "red" , linewidth = 0.9) +
  facet_wrap(.~ Model, scales = "free_x", nrow = 1) +
  #geom_ribbon(aes(ymin = lb_bands, ymax = ub_bands), color = "red" , linewidth = 0.8, alpha = 0.1) +
  geom_point( data = accum_df, aes(x = x, y = n_feat), color="black", shape = 1, size = 0.5) +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1)
ggsave(filename = "R_script_paper/Paper_plots/rarefaction_dune_eb_MM.pdf", width = 8, height = 3.5, dpi = 300, units = "in", device='pdf')

# 3 variances 
n_rare <- n
eb_fit_NegBinBB_EFPF_1_rare <- eb_NegBinBB_list_variances[[1]]
eb_fit_NegBinBB_EFPF_2_rare <- eb_NegBinBB_list_variances[[2]]
eb_fit_NegBinBB_EFPF_3_rare <- eb_NegBinBB_list_variances[[3]]

eb_fit_NegBinBB_MM_1_rare <- eb_NegBinBB_MM_list_variances[[1]]
eb_fit_NegBinBB_MM_2_rare <- eb_NegBinBB_MM_list_variances[[2]]
eb_fit_NegBinBB_MM_3_rare <- eb_NegBinBB_MM_list_variances[[3]]

accum_df <- tibble( x = 0:n_rare,
                    n_feat = c(0,rarefaction(data_mat[1:n_rare,], n_reorderings = 10)))



rare_NegBinBB_EFPF_1 <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_fit_NegBinBB_EFPF_1_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = paste0("EFPF:",cfrs[1]," x"),
             x = c(1:n_rare,0))

rare_NegBinBB_MM_1 <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_fit_NegBinBB_MM_1_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = paste0("MM:",cfrs[1]," x"),
             x = c(1:n_rare,0))

rare_NegBinBB_EFPF_2 <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_fit_NegBinBB_EFPF_2_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = paste0("EFPF:", cfrs[2]," x"),
             x = c(1:n_rare,0))

rare_NegBinBB_MM_2 <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_fit_NegBinBB_MM_2_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = paste0("MM:", cfrs[2]," x"),
             x = c(1:n_rare,0))

rare_NegBinBB_EFPF_3 <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_fit_NegBinBB_EFPF_3_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = paste0("EFPF:", cfrs[3]," x"),
             x = c(1:n_rare,0))

rare_NegBinBB_MM_3 <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_fit_NegBinBB_MM_3_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = paste0("MM:", cfrs[3]," x"),
             x = c(1:n_rare,0))

df_rare <- rbind(#rare_NegBinBB_EFPF_1, rare_NegBinBB_EFPF_2, rare_NegBinBB_EFPF_3,
                 rare_NegBinBB_MM_1, rare_NegBinBB_MM_2, rare_NegBinBB_MM_3)

df_rare$Model <- factor(df_rare$Model)

ggplot(df_rare, aes(x = x, y = means, color = Model)) +
  geom_line(linetype = "solid", color = "red" , linewidth = 0.9) +
  facet_wrap(.~ Model, scales = "free_x", nrow = 1) +
  #geom_ribbon(aes(ymin = lb_bands, ymax = ub_bands), color = "red" , linewidth = 0.8, alpha = 0.1) +
  geom_point( data = accum_df, aes(x = x, y = n_feat), color="black", shape = 1, size = 0.5) +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1)
ggsave(filename = "R_script_paper/Paper_plots/rarefaction_variances_dune_eb_MM.pdf", width = 8, height = 3.5, dpi = 300, units = "in", device='pdf')


###### Check on K_n_r -----
n_knr <- n
# eb_fit_NegBinBB_EFPF_1_knr <- eb_NegBinBB_list_variances[[1]]
# eb_fit_NegBinBB_EFPF_2_knr <- eb_NegBinBB_list_variances[[2]]
# eb_fit_NegBinBB_EFPF_3_knr <- eb_NegBinBB_list_variances[[3]]
# eb_fit_NegBinBB_EFPF_4_knr <- eb_NegBinBB_list_variances[[4]]
# 
# eb_fit_NegBinBB_MM_1_knr <- eb_NegBinBB_MM_list_variances[[1]]
# eb_fit_NegBinBB_MM_2_knr <- eb_NegBinBB_MM_list_variances[[2]]
# eb_fit_NegBinBB_MM_3_knr <- eb_NegBinBB_MM_list_variances[[3]]
# eb_fit_NegBinBB_MM_4_knr <- eb_NegBinBB_MM_list_variances[[4]]

observed_K_n_r <- tibble( r = 1:n_knr,
                    k_n_r = K_n_r(data_mat[1:n_knr,], n_reorderings = 1)[[paste0('N = ', n_knr)]])

K_n_r_PoissonBB_EFPF <- tibble( lambda_est = unname(unlist(
  K_n_r(eb_PoissonBB_fit, n =n_knr)[[paste0('N = ', n_knr)]]$lambda_est ))) %>%
  rename(means = lambda_est) %>%
  add_column(Model = "PoissonBB_EFPF",
             r = 1:n_knr)

K_n_r_PoissonBB_MM_biased <- tibble( lambda_est = unname(unlist(
  K_n_r(eb_PoissonBB_MM_biased_fit, n =n_knr)[[paste0('N = ', n_knr)]]$lambda_est ))) %>%
  rename(means = lambda_est) %>%
  add_column(Model = "PoissonBB_MM_biased",
             r = 1:n_knr)

K_n_r_PoissonBB_MM_cens <- tibble( lambda_est = unname(unlist(
  K_n_r(eb_PoissonBB_MM_cens_fit, n =n_knr)[[paste0('N = ', n_knr)]]$lambda_est ))) %>%
  rename(means = lambda_est) %>%
  add_column(Model = "PoissonBB_MM_cens",
             r = 1:n_knr)

# K_n_r_NegBinBB_EFPF_1 <- tibble( mu0_est = unname(unlist(
#   K_n_r(object = eb_fit_NegBinBB_EFPF_1_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
#   rename(means = mu0_est) %>%
#   add_column(Model = paste0("EFPF, ",cfrs[1], "x"),
#              r = 1:n_knr)
# 
# K_n_r_NegBinBB_MM_1 <- tibble( mu0_est = unname(unlist(
#   K_n_r(object = eb_fit_NegBinBB_MM_1_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
#   rename(means = mu0_est) %>%
#   add_column(Model = paste0("MM, ",cfrs[1], "x"),
#              r = 1:n_knr)
# 
# K_n_r_NegBinBB_EFPF_2 <- tibble( mu0_est = unname(unlist(
#   K_n_r(object = eb_fit_NegBinBB_EFPF_2_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
#   rename(means = mu0_est) %>%
#   add_column(Model = paste0("EFPF, ",cfrs[2], "x"),
#              r = 1:n_knr)
# 
# K_n_r_NegBinBB_MM_2 <- tibble( mu0_est = unname(unlist(
#   K_n_r(object = eb_fit_NegBinBB_MM_2_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
#   rename(means = mu0_est) %>%
#   add_column(Model = paste0("MM, ",cfrs[2], "x"),
#              r = 1:n_knr)
# 
# K_n_r_NegBinBB_EFPF_3 <- tibble( mu0_est = unname(unlist(
#   K_n_r(object = eb_fit_NegBinBB_EFPF_3_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
#   rename(means = mu0_est) %>%
#   add_column(Model = paste0("EFPF, ",cfrs[3], "x"),
#              r = 1:n_knr)
# 
# K_n_r_NegBinBB_MM_3 <- tibble( mu0_est = unname(unlist(
#   K_n_r(object = eb_fit_NegBinBB_MM_3_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
#   rename(means = mu0_est) %>%
#   add_column(Model = paste0("MM, ",cfrs[3], "x"),
#              r = 1:n_knr)
# 
# 
# K_n_r_NegBinBB_EFPF_4 <- tibble( mu0_est = unname(unlist(
#   K_n_r(object = eb_fit_NegBinBB_EFPF_4_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
#   rename(means = mu0_est) %>%
#   add_column(Model = paste0("EFPF, ",cfrs[4], "x"),
#              r = 1:n_knr)
# 
# K_n_r_NegBinBB_MM_4 <- tibble( mu0_est = unname(unlist(
#   K_n_r(object = eb_fit_NegBinBB_MM_4_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
#   rename(means = mu0_est) %>%
#   add_column(Model = paste0("MM, ",cfrs[4], "x"),
#              r = 1:n_knr)
# 
K_n_r_GammaIBP_EFPF <- tibble( mu0_est = unname(unlist(
  K_n_r(eb_GammaIBP_EFPF_fit, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
  rename(means = mu0_est) %>%
  add_column(Model = "GammaIBP_EFPF",
             r = 1:n_knr)

K_n_r_GammaIBP_MM <- tibble( mu0_est = unname(unlist(
  K_n_r(eb_GammaIBP_MM_fit, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
  rename(means = mu0_est) %>%
  add_column(Model = "GammaIBP_MM",
             r = 1:n_knr)

df_K_n_r <- rbind(K_n_r_PoissonBB_EFPF, K_n_r_PoissonBB_MM_biased, K_n_r_PoissonBB_MM_cens,
                  #K_n_r_NegBinBB_EFPF_1, K_n_r_NegBinBB_EFPF_2, K_n_r_NegBinBB_EFPF_3, K_n_r_NegBinBB_EFPF_4,
                  # K_n_r_NegBinBB_MM_1, K_n_r_NegBinBB_MM_2, K_n_r_NegBinBB_MM_3, K_n_r_NegBinBB_MM_4)
                  K_n_r_GammaIBP_EFPF, K_n_r_GammaIBP_MM)

                  
df_K_n_r$Model <- factor(df_K_n_r$Model)

r_positive <- observed_K_n_r %>%
  filter(k_n_r > 0) %>%
  select(r)

df_K_n_r_plot <- df_K_n_r %>%
  filter(r %in% c(r_positive$r))

observed_K_n_r_plot <- observed_K_n_r %>%
  filter(r %in% c(r_positive$r))

ggplot(df_K_n_r_plot, aes(x = r, y = means, color = Model)) +
  geom_line( linetype = "dashed") +
  geom_point( data = observed_K_n_r_plot, aes(x = r, y = k_n_r), color="black", size = 1) +
  scale_y_log10() +
  xlab("r") + ylab(expression(K[paste("n,r")])) + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1)
ggsave(filename = "R_script_paper/Paper_plots/knr_dune_eb_MM.pdf", width = 5, height = 4.5, dpi = 300, units = "in", device='pdf')





# 3) Estimate the total richness ----

lambda_prime <- total_richness(eb_PoissonBB_fit)$lambda_post
lb_PoissonBB <- qpois(0.005, lambda_prime, lower.tail = TRUE, log.p = FALSE)
ub_PoissonBB <- qpois(0.995, lambda_prime, lower.tail = TRUE, log.p = FALSE)

mu0_prime <- total_richness(eb_NegBinBB_fit)$mu0_post
n0_prime <- total_richness(eb_NegBinBB_fit)$n0_post
p_prime = 1/(mu0_prime/n0_prime + 1)
lb_NegBinBB <- qnbinom(0.005, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE)
ub_NegBinBB <- qnbinom(0.995, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE) 


bounds <- tibble(lb = min(lb_PoissonBB, lb_NegBinBB) + Kn, 
            ub = max(ub_PoissonBB, ub_NegBinBB) + Kn )

dens_richness_PoissonBB <- tibble( x = bounds$lb: bounds$ub) %>%
  mutate( y = dpois(x - Kn , lambda = lambda_prime)) %>%
  add_column(Model = "PoissonBB")

dens_richness_NegBinBB <- tibble( x = bounds$lb : bounds$ub) %>%
  mutate( y = dnbinom(x - Kn, size = n0_prime, 
                      prob = p_prime)) %>%
  add_column(Model = "NegBinBB")

dens_richnesses <- rbind(dens_richness_PoissonBB,
                         dens_richness_NegBinBB)


ggplot(dens_richnesses, aes(x = x, y = y, color = Model)) +
  geom_line() +
  theme_light() +
  facet_wrap(~"Richness") +
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  xlab("# distinct features") + rremove("ylab") +
  scale_color_tableau() +
  theme(aspect.ratio = 1)
ggsave(filename = "R_script_paper/Paper_plots/richness_dune_eb.pdf", width = 3.8, height = 3.8, dpi = 300, units = "in", device='pdf')


# 3) Plot Extrapolation - EB version

# Extract accumulation curve of the observed sample (or average accumulation)
M = 200

accum_df <- tibble( x = 0:n,
                    n_feat = c(0,rarefaction(data_mat, n_reorderings = 20)))


extr_PoissonBB_df <- tibble(lambda = unname(unlist( 
  extrapolation(object = eb_PoissonBB_fit, M = M, seed = seed)$lambda_post)),
  Kn = rep(Kn, each = M)) %>%
  mutate(lb = qpois(0.025, lambda, lower.tail = TRUE, log.p = FALSE),
         ub = qpois(0.975, lambda, lower.tail = TRUE, log.p = FALSE)) %>%
  rename(means = lambda) %>%
  add_row(means = 0, lb = 0, ub = 0, Kn = Kn) %>%
  mutate(means = means + Kn, lb = lb + Kn, ub = ub + Kn) %>%
  add_column(x = c((n+1):(n+M), n),
             Model = "PoissonBB")

extr_PoissonBB_df$x <- as.integer(extr_PoissonBB_df$x)
extr_PoissonBB_df <- extr_PoissonBB_df %>%
  select(means, lb, ub, x, Model)

extr_NegBinBB_df <- tibble(mu0 = unname(unlist(
  extrapolation(object = eb_NegBinBB_fit, M = M, seed = seed)$mu0_post )),
  n0 = unname(unlist( extrapolation(object = eb_NegBinBB_fit, M = M, seed = seed)$n0_post )),
  Kn = rep(Kn, each = M)) %>%
  mutate(p = 1/(mu0/n0 + 1),
         lb = qnbinom(0.025, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE),
         ub = qnbinom(0.975, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE)) %>%
  rename(means = mu0) %>%
  add_row(means = 0, lb = 0, ub = 0, Kn = Kn) %>%
  mutate(means = means + Kn, lb = lb + Kn, ub = ub + Kn) %>%
  add_column(x = c((n+1):(M+n), n),
             Model = "NegBinBB")

extr_NegBinBB_df$x <- as.integer(extr_NegBinBB_df$x)
extr_NegBinBB_df <- extr_NegBinBB_df %>%
  select(means, lb, ub, x, Model)

extr_GammaIBP_df <- tibble(mu0 = unname(unlist(
  extrapolation(object = eb_GammaIBP_fit, M = M, seed = seed)$mu0_post )),
  n0 = unname(unlist( extrapolation(object = eb_GammaIBP_fit, M = M, seed = seed)$n0_post )),
  Kn = rep(Kn, each = M)) %>%
  mutate(p = 1/(mu0/n0 + 1),
         lb = qnbinom(0.025, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE),
         ub = qnbinom(0.975, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE)) %>%
  rename(means = mu0) %>%
  add_row(means = 0, lb = 0, ub = 0, Kn = Kn) %>%
  mutate(means = means + Kn, lb = lb + Kn, ub = ub + Kn) %>%
  add_column(x = c((n+1):(M+n), n),
             Model = "GammaIBP")

extr_GammaIBP_df$x <- as.integer(extr_GammaIBP_df$x)
extr_GammaIBP_df <- extr_GammaIBP_df %>%
  select(means, lb, ub, x, Model)


extr_all_df <- rbind(extr_PoissonBB_df, extr_NegBinBB_df, extr_GammaIBP_df)

ggplot(extr_all_df, aes(x, means, color = Model)) +
  geom_line(linetype = "dashed", linewidth = 1) +
  geom_point( data = accum_df, aes(x = x, y = n_feat),
              color="black", shape = 1, size = 1) +
  #geom_line(data = df_extr_GT_long, aes(t, value)) +
  #geom_line(data = df_extr_Chao_long, aes(t, medians), linetype = "dashed", linewidth = 1) +
  #facet_wrap(. ~ n_train,labeller = labeller(n_train = n_train.labs ),  scales = "free_x") +
  geom_vline(aes(xintercept = n) , linetype = "dashed", color = "grey") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  facet_wrap(~"Extrapolation") +
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1)
ggsave(filename = "R_script_paper/Paper_plots/extr_dune_eb.pdf", width = 3.8, height = 3.8, dpi = 300, units = "in", device='pdf')



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








