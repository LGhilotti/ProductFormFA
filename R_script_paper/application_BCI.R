rm(list=ls())

library(ProductFormFA)
library(tidyverse)
library(scales)
library(ggpubr)
library(ggthemes)

source("R_script_paper/Routine_Chao.R")
source("R_script_paper/utils.R")

# Fit models -----

#data(birds, package="jSDM")
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
  geom_point(color="black", shape = 21, size = 1) + 
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) + 
  scale_color_tableau()
ggsave(filename = "R_script_paper/Paper_plots/accumulation_BCI.pdf", width = 4, height = 4, dpi = 300, units = "in", device='pdf')


# 2) Run the models on the data

vars_fct_NegBinBB <- c(10, 1000)
vars_GammaIBP <- c(1, 1000)

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



# Model Checking ----

# 0.B) Check on rarefaction
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

# for plot
rare_EFPF_mixtureBB <- rare_EFPF_NegBinBB %>%
  mutate(Model = "PoissonBB/NegBinBB")

df_rare <- rbind(rare_EFPF_mixtureBB,
                 rare_EFPF_GammaIBP)

df_rare$Model <- factor(df_rare$Model,
                        levels = c("PoissonBB/NegBinBB", "GammaIBP"))

ggplot(df_rare, aes(x = x, y = means, color = Model)) +
  geom_line(linetype = "dashed", linewidth = 1) + 
  #facet_wrap(.~ Model, scales = "free_x", nrow = 1) +
  #geom_ribbon(aes(ymin = lb_bands, ymax = ub_bands), color = "red" , linewidth = 0.8, alpha = 0.1) +
  geom_point( data = accum_df, aes(x = x, y = n_feat), color="black", shape = 21, size = 0.5) +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) + 
  scale_color_tableau()
ggsave(filename = "R_script_paper/Paper_plots/rarefaction_BCI_eb_EFPF.pdf", width = 4, height = 4, dpi = 300, units = "in", device='pdf')


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
  filter(r < 15)

df_K_n_r_plot <- df_K_n_r %>%
  filter(r %in% c(r_positive$r))

observed_K_n_r_plot <- observed_K_n_r %>%
  filter(r %in% c(r_positive$r))

ggplot(df_K_n_r_plot, aes(x = r, y = means, color = Model)) +
  geom_line( linetype = "dashed") +
  geom_point( data = observed_K_n_r_plot, aes(x = r, y = k_n_r), color="black", shape = 21, size = 1) +
  scale_y_log10() +
  xlab("r") + ylab(expression(m[r])) + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) + 
  scale_color_tableau()
ggsave(filename = "R_script_paper/Paper_plots/knr_BCI_eb_EFPF.pdf", width = 4, height = 4, dpi = 300, units = "in", device='pdf')



# Prediction ----

## Richness estimation -----

# PoissonBB
params_richness_EFPF_PoissonBB <- tibble( lambda_prime = 
                                            total_richness(eb_EFPF_fit_PoissonBB)$lambda_post) %>%
  add_column(Model = "PoissonBB") %>%
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
    add_column(Model = paste0("NegBinBB x", var_fct_NegBinBB)) %>%
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
  add_column(Model = "PoissonBB")

# NegBinBB
dens_richness_NegBinBB <- tibble(x = integer(), y = numeric(),
                                 Model = character() )

for (var_fct_NegBinBB in vars_fct_NegBinBB){
  
  params_richness_EFPF_NegBinBB_var <- params_richness_EFPF_NegBinBB %>%
    filter(Model == paste0("NegBinBB x", var_fct_NegBinBB))
  
  dens_richness_NegBinBB_var <- tibble( x = bounds$lb : bounds$ub) %>%
    mutate( y = dnbinom(x - Kn, size = params_richness_EFPF_NegBinBB_var$n0_prime, 
                        prob = params_richness_EFPF_NegBinBB_var$p_prime)) %>%
    add_column(Model = paste0("NegBinBB x", var_fct_NegBinBB))
  
  dens_richness_NegBinBB <- bind_rows(dens_richness_NegBinBB,
                                      dens_richness_NegBinBB_var)
  
}



dens_richnesses <- rbind(dens_richness_PoissonBB,
                         dens_richness_NegBinBB) 

# dens_richnesses$Model[dens_richnesses$Model == paste0("NegBinBB x", vars_fct_NegBinBB)] <-
#   "NegBinBB"

dens_richnesses$Model <- factor(dens_richnesses$Model, 
                                levels = c("PoissonBB",# "NegBinBB"))
                                           paste0("NegBinBB x", vars_fct_NegBinBB)))



ggplot(dens_richnesses, aes(x = x, y = y, color = Model)) +
  #geom_col(position = "identity", alpha = 0.1) +
  geom_line() +
  theme_light() +
  #facet_wrap(~"Richness") +
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  xlab("# distinct features") + rremove("ylab") +
  scale_color_tableau() +
  theme(aspect.ratio = 1)
ggsave(filename = "R_script_paper/Paper_plots/richness_BCI_eb_EFPF.pdf", width = 4.5, height = 4.5, dpi = 300, units = "in", device='pdf')


# mean and variace species richness

# PoissonBB
rich_pars <- params_richness_EFPF_PoissonBB 

print(paste0("mean richness = ", rich_pars$lambda_prime + Kn ))
Kn + qpois(c(0.025, 0.975), lambda = rich_pars$lambda_prime)


# NegBinBB
rich_pars <- params_richness_EFPF_NegBinBB %>%
  filter(Model == "NegBinBB x1000")

print(paste0("mean richness = ", rich_pars$mu0_prime + Kn ))
Kn + qnbinom(c(0.025, 0.975), size = rich_pars$n0_prime, prob = rich_pars$p_prime)


## Extrapolation (EFPF version) -----

# Extract accumulation curve of the observed sample (or average accumulation)
M = 100

accum_df <- tibble( x = 0:n,
                    n_feat = c(0,rarefaction(data_mat, n_reorderings = 20)))


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
             Model = "PoissonBB")

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
               Model = paste0("NegBinBB x", var_fct_NegBinBB))
  
  extr_EFPF_NegBinBB_df_var$x <- as.integer(extr_EFPF_NegBinBB_df_var$x)
  extr_EFPF_NegBinBB_df_var <- extr_EFPF_NegBinBB_df_var %>%
    select(means, lb, ub, x, Model)
  
  extr_EFPF_NegBinBB_df <- bind_rows(extr_EFPF_NegBinBB_df, 
                                     extr_EFPF_NegBinBB_df_var)
}

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
               Model = paste0("GammaIBP, var = ", var_GammaIBP))
  
  extr_GammaIBP_df_var$x <- as.integer(extr_GammaIBP_df_var$x)
  extr_GammaIBP_df_var <- extr_GammaIBP_df_var %>%
    select(means, lb, ub, x, Model)
  
  extr_EFPF_GammaIBP_df <- bind_rows(extr_EFPF_GammaIBP_df, 
                                     extr_GammaIBP_df_var)
}

extr_EFPF_PoissonBB_df <- extr_EFPF_PoissonBB_df %>%
  add_column(Model_gen = "PoissonBB/NegBinBB")
extr_EFPF_NegBinBB_df <- extr_EFPF_NegBinBB_df %>%
  add_column(Model_gen = "PoissonBB/NegBinBB") 
#extr_EFPF_NegBinBB_df$Model <- "NegBinBB"
extr_EFPF_GammaIBP_df <- extr_EFPF_GammaIBP_df %>%
  add_column(Model_gen = "GammaIBP")

extr_all_df <- rbind(extr_EFPF_PoissonBB_df, 
                     extr_EFPF_NegBinBB_df)#,
                     #extr_EFPF_GammaIBP_df)


extr_all_df$Model <- factor(extr_all_df$Model,
                            levels = c("PoissonBB",#"NegBinBB"))
                                       paste0("NegBinBB x", vars_fct_NegBinBB)))

extr_all_df$Model_gen <- factor(extr_all_df$Model_gen,
                                levels = c("PoissonBB/NegBinBB", "GammaIBP"))

ggplot(extr_all_df, aes(x, means, color = Model)) +
  geom_line(linetype = "dashed", linewidth = 1) +
  #facet_wrap(. ~ Model_gen,  scales = "free_x") +
  geom_point( data = accum_df, aes(x = x, y = n_feat),
              color="black", shape = 1, size = 1) +
  geom_ribbon(aes(ymin = lb, ymax = ub, color = Model), linewidth = 0.8, alpha = 0) +
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
ggsave(filename = "R_script_paper/Paper_plots/extr_BCI_eb_EFPF.pdf", width = 4.5, height = 4.5, dpi = 300, units = "in", device='pdf')


extr_EFPF_NegBinBB_df %>%
  filter(Model == "NegBinBB x1000",
         x %in% c(n + 1, n + 10, n + 100)) %>%
  mutate(means_new = means - Kn,
         lb_new = lb -Kn,
         ub_new = ub - Kn)

extr_EFPF_PoissonBB_df %>%
  filter(x %in% c(n + 1, n + 10, n + 100)) %>%
  mutate(means_new = means - Kn,
         lb_new = lb -Kn,
         ub_new = ub - Kn)

## Extrapolation (prior approach)- to be done ----

PoissonBB_extr <- extrapolation(object = PoissonBB_fit, M = M) 
NegBinBB_extr <- extrapolation(object = NegBinBB_fit, M = M) 


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








