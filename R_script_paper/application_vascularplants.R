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
vars_GammaIBP <- c(100, 1000, 10000) # c(0.01, 100) - values in the first manuscript

# Initial parameters for optimization
eb_init_BB <- list(alpha = -10, s = 100, Nhat_prime = 200)
eb_known_BB <- list()

eb_init_IBP <- list(alpha = 0.4, s = 2, Gamma = 10)
eb_known_IBP <- list()

eb_params_obj_BB <- eb_params(model = "BB", 
                              init = eb_init_BB, known = eb_known_BB )
eb_params_obj_IBP <- eb_params(model = "IBP", 
                               init = eb_init_IBP, known = eb_known_IBP )

# classicBB
eb_EFPF_fit_classicBB <- GibbsFA_eb(feature_matrix = data_mat, 
                                    model = "classicBB_eb", 
                                    type = "EFPF",
                                    eb_params =  eb_params_obj_BB)

# PoissonBB
eb_EFPF_fit_PoissonBB <- GibbsFA_eb(feature_matrix = data_mat, 
                                    model = "PoissonBB_eb", 
                                    type = "EFPF",
                                    eb_params =  eb_params_obj_BB)

# NegBinBB
list_eb_EFPF_fit_NegBinBB <-  vector(mode = "list", length = length(vars_fct_NegBinBB))
names(list_eb_EFPF_fit_NegBinBB) <- paste0("var_fct.", vars_fct_NegBinBB)

for (var_fct_NegBinBB in vars_fct_NegBinBB){
  
  list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]] <- 
    GibbsFA_eb(feature_matrix = data_mat,
               model = "NegBinBB_eb", type = "EFPF",
               eb_params =  eb_params_obj_BB, 
               var_fct = var_fct_NegBinBB)
  
}

# classicIBP
eb_EFPF_fit_classicIBP <- GibbsFA_eb(feature_matrix = data_mat, 
                                     model = "classicIBP_eb", 
                                     type = "EFPF",
                                     eb_params =  eb_params_obj_IBP)

# GammaIBP
list_eb_EFPF_fit_GammaIBP <-  vector(mode = "list", length = length(vars_GammaIBP))
names(list_eb_EFPF_fit_GammaIBP) <- paste0("var.", vars_GammaIBP)

for (var_GammaIBP in vars_GammaIBP){
  
  list_eb_EFPF_fit_GammaIBP[[paste0("var.", var_GammaIBP)]] <-
    GibbsFA_eb(feature_matrix = data_mat,
               model = "GammaIBP_eb", type = "EFPF",
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
  filter(r < 50)

df_K_n_r_plot <- df_K_n_r %>%
  filter(r %in% c(r_positive$r))

observed_K_n_r_plot <- observed_K_n_r %>%
  filter(r %in% c(r_positive$r))

# for plot
ggplot(observed_K_n_r_plot,  aes(x = r, y = k_n_r)) +
  geom_point(color="black", shape = 19, size = 1) +
  geom_line( data = df_K_n_r_plot, aes(x = r, y = means, color = Model), linetype = "dashed") +
  scale_y_log10() +
  scale_x_log10() +
  xlab("r") + ylab(expression(m[r])) + 
  theme_light() + 
  theme(legend.position = "top") +
  #scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) + 
  scale_color_tableau(
    labels = c(
      "PoissonBB/NegBinBB" = "Mixtures of BBs",
      "GammaIBP" = "Mixtures of IBPs"
    )
  )
#ggsave(filename = "R_script_paper/Paper_plots/knr_Plants_eb_EFPF.pdf", width = 4, height = 4, dpi = 300, units = "in", device='pdf')


## Formal model-checking via AIC/BIC -------
AICs_list <- vector("list", length = 0)

AICs_list[["classicBB"]] <- compute_AICs_BICs(eb_EFPF_fit_classicBB)$AIC
AICs_list[["PoissonBB"]] <- compute_AICs_BICs(eb_EFPF_fit_PoissonBB)$AIC
for (var_fct_NegBinBB in vars_fct_NegBinBB){
  AICs_list[[paste0("NegBinBB.var_fct.", var_fct_NegBinBB)]] <- compute_AICs_BICs(
    list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]]
  )$AIC
}
AICs_list[["classicIBP"]] <- compute_AICs_BICs(eb_EFPF_fit_classicIBP)$AIC
for (var_GammaIBP in vars_GammaIBP){
  AICs_list[[paste0("GammaIBP.var.", var_GammaIBP)]] <- compute_AICs_BICs(
    list_eb_EFPF_fit_GammaIBP[[paste0("var.", var_GammaIBP)]]
  )$AIC
}

print(AICs_list) 
# PoissonBB/NegBinBB are worse than GammaIBP (higher AIC)

# Max-log_efpf
max_log_efpf_list <- vector("list", length = 0)

max_log_efpf_list[["classicBB"]] <- compute_AICs_BICs(eb_EFPF_fit_classicBB)$max_log_efpf
max_log_efpf_list[["PoissonBB"]] <- compute_AICs_BICs(eb_EFPF_fit_PoissonBB)$max_log_efpf
for (var_fct_NegBinBB in vars_fct_NegBinBB){
  max_log_efpf_list[[paste0("NegBinBB.var_fct.", var_fct_NegBinBB)]] <- compute_AICs_BICs(
    list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]]
  )$max_log_efpf
}
max_log_efpf_list[["classicIBP"]] <- compute_AICs_BICs(eb_EFPF_fit_classicIBP)$max_log_efpf
for (var_GammaIBP in vars_GammaIBP){
  max_log_efpf_list[[paste0("GammaIBP.var.", var_GammaIBP)]] <- compute_AICs_BICs(
    list_eb_EFPF_fit_GammaIBP[[paste0("var.", var_GammaIBP)]]
  )$max_log_efpf
}

max_log_efpf_df <- data.frame(
  Model = names(max_log_efpf_list),
  Max_log_epfp = unlist(max_log_efpf_list),
  row.names = NULL
)
print(max_log_efpf_df)
#write.csv(max_log_efpf_df, file = "max_log_efpf_df.csv", row.names = FALSE)



## Rarefaction and Knr plots with credible bands for best class of mixtures -----------

### Rarefaction intervals for GammaIBP -----
n_rare <- n

accum_df <- tibble( x = 0:n_rare,
                    n_feat = c(0,rarefaction(data_mat[1:n_rare,], n_reorderings = 20)))




# GammaIBP
rare_EFPF_GammaIBP_df <- tibble(means = numeric(), 
                                lb = numeric(), ub = numeric(),
                                x = integer(), Model = character())

for (var_GammaIBP in vars_GammaIBP){
  
  eb_EFPF_fit_GammaIBP_var <- 
    list_eb_EFPF_fit_GammaIBP[[paste0("var.", var_GammaIBP)]]
  
  rare_EFPF_GammaIBP_df_var <- tibble( mu0_post = unname(unlist(
    rarefaction(object = eb_EFPF_fit_GammaIBP_var, seed = seed)$mu0_post )),
    n0_post = unname(unlist(
      rarefaction(object = eb_EFPF_fit_GammaIBP_var, seed = seed)$n0_post ))) %>%
    mutate(p_post = 1/(mu0_post/n0_post + 1),
           lb = qnbinom(0.025, size = n0_post, prob = p_post, lower.tail = TRUE, log.p = FALSE),
           ub = qnbinom(0.975, size = n0_post, prob = p_post, lower.tail = TRUE, log.p = FALSE)) %>%
    rename(means = mu0_post) %>%
    add_row(means = 0, lb = 0, ub = 0) %>%
    add_column(Model = paste0("Gamma IBP, Variance: ", var_GammaIBP),
               x = c(1:n_rare,0))
  
  rare_EFPF_GammaIBP_df_var$x <- as.integer(rare_EFPF_GammaIBP_df_var$x)
  rare_EFPF_GammaIBP_df_var <- rare_EFPF_GammaIBP_df_var %>%
    select(means, lb, ub, x, Model)
  
  rare_EFPF_GammaIBP_df <- bind_rows(rare_EFPF_GammaIBP_df, 
                                     rare_EFPF_GammaIBP_df_var)
  
}


rare_all_df <- rare_EFPF_GammaIBP_df

rare_all_df$Model <- factor(rare_all_df$Model, 
                            levels = paste0("Gamma IBP, Variance: ", vars_GammaIBP))


# for plot
plot_ribbons_rare <- ggplot() +
  geom_point(data = accum_df, aes(x = x, y = n_feat),
             color="black", shape = 18, size = 1) +
  geom_ribbon(data = rare_all_df, aes(x = x, ymin = lb, ymax = ub, fill = Model), color = NA, alpha = 0.4) +
  scale_fill_manual(values = c(
    "Gamma IBP, Variance: 100" = "grey10",
    "Gamma IBP, Variance: 1000" = "grey50",
    "Gamma IBP, Variance: 10000" = "grey80")) +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) +
  scale_color_tableau()

plot_ribbons_rare


### Knr intervals for GammaIBP --------
n_knr <- n

observed_K_n_r <- tibble( r = 1:n_knr,
                          k_n_r = K_n_r(data_mat[1:n_knr,], n_reorderings = 1)[[paste0('N = ', n_knr)]])


# GammaIBP
knr_EFPF_GammaIBP_df <- tibble(means = numeric(), 
                               lb = numeric(), ub = numeric(),
                               r = integer(), Model = character())

for (var_GammaIBP in vars_GammaIBP){
  
  eb_EFPF_fit_GammaIBP_var <- 
    list_eb_EFPF_fit_GammaIBP[[paste0("var.", var_GammaIBP)]]
  
  knr_EFPF_GammaIBP_df_var <- tibble( mu0_est = unname(unlist(
    K_n_r(object = eb_EFPF_fit_GammaIBP_var, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est )),
    n0_est = unname(unlist(
      K_n_r(object = eb_EFPF_fit_GammaIBP_var, n = n_knr)[[paste0('N = ', n_knr)]]$n0_est ))) %>%
    mutate(p_est = 1/(mu0_est/n0_est + 1),
           lb = qnbinom(0.025, size = n0_est, prob = p_est, lower.tail = TRUE, log.p = FALSE),
           ub = qnbinom(0.975, size = n0_est, prob = p_est, lower.tail = TRUE, log.p = FALSE)) %>%
    rename(means = mu0_est) %>%
    add_column(Model = paste0("Gamma IBP, Variance: ", var_GammaIBP),
               r = 1:n_knr)
  
  knr_EFPF_GammaIBP_df_var$r <- as.integer(knr_EFPF_GammaIBP_df_var$r)
  knr_EFPF_GammaIBP_df_var <- knr_EFPF_GammaIBP_df_var %>%
    select(means, lb, ub, r, Model)
  
  knr_EFPF_GammaIBP_df <- bind_rows(knr_EFPF_GammaIBP_df, 
                                    knr_EFPF_GammaIBP_df_var)
  
}


knr_all_df <- knr_EFPF_GammaIBP_df

knr_all_df$Model <- factor(knr_all_df$Model,
                           levels = paste0("Gamma IBP, Variance: ", vars_GammaIBP))

r_positive <- observed_K_n_r %>%
  filter(k_n_r > 0) %>%
  select(r) %>%
  filter(r < 50)

knr_all_df_plot <- knr_all_df %>%
  filter(r %in% c(r_positive$r)) %>%
  mutate(lb = ifelse(lb == 0, 8e-1, lb))
# %>% filter(Model %in% c("Poisson BB", "NegBinomial BB x10"))

observed_K_n_r_plot <- observed_K_n_r %>%
  filter(r %in% c(r_positive$r))


# for plot
plot_ribbons_knr <- ggplot() +
  geom_point(data = observed_K_n_r_plot, aes(x = r, y = k_n_r),
             color="black", shape = 19, size = 1.5) +
  geom_ribbon(data = knr_all_df_plot, aes(x = r, ymin = lb, ymax = ub, fill = Model), color = NA, alpha = 0.4) +
  scale_fill_manual(values = c(
    "Gamma IBP, Variance: 100" = "grey10",
    "Gamma IBP, Variance: 1000" = "grey50",
    "Gamma IBP, Variance: 10000" = "grey80")) +
  scale_y_log10() +
  scale_x_log10() +
  xlab("r") + ylab(expression(m[r])) + 
  theme_light() + 
  theme(legend.position = "top") +
  #scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) +
  scale_color_tableau()


plot_ribbons_knr





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
if (!file.exists("R_script_paper/fullybayes_Plants_fit_estimate_GammaIBP.RData")) {
  
  vars_GammaIBP_bayes <- c(100, 1000, 10000) # c(0.01, 100) - values in the first manuscript
  
  list_prior_fit_GammaIBP <-  vector(mode = "list", length = length(vars_GammaIBP_bayes))
  names(list_prior_fit_GammaIBP) <- paste0("var.", vars_GammaIBP_bayes)
  
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
  for (var_GammaIBP in vars_GammaIBP_bayes){
    
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
  save(list_prior_fit_GammaIBP, vars_GammaIBP_bayes, file =  "R_script_paper/fullybayes_Plants_fit_estimate_GammaIBP.RData")
  
}

# We run also for NegBinBB + prior (in order to check with BF other than visual check)
if (!file.exists("R_script_paper/fullybayes_Plants_fit_estimate_NegBinBBcompetitor.RData")) {
  
  vars_fct_NegBinBB_bayes <- c(10, 1000)
  
  list_prior_fit_NegBinBB <-  vector(mode = "list", length = length(vars_fct_NegBinBB_bayes))
  names(list_prior_fit_NegBinBB) <- paste0("var_fct.", vars_fct_NegBinBB_bayes)
  
  # Initialization and MCMC setting 
  mcmcparams_NegBinBB <- list(tau = 0.1, 
                              S = 5*10^4, n_burnin = 5*10^3, thin = 2)
  mcmcparams_obj_NegBinBB <- mcmcparameters(model = "NegBinBB", mcmcparams = mcmcparams_NegBinBB)
  
  init_NegBinBB <- list(alpha_0 = - 0.01, s_0 = 3)
  init_obj_NegBinBB <- initialization(model = "NegBinBB", init = init_NegBinBB )
  
  # EB estimates
  small_val <- 10^(-4) # 10^(-3)
  alpha_eb <- list_eb_EFPF_fit_NegBinBB[[1]]$alpha
  theta_eb <- list_eb_EFPF_fit_NegBinBB[[1]]$theta
  
  s_eb <- alpha_eb + theta_eb
  
  print(paste0("Prior variance of -alpha: ", - alpha_eb/small_val  ))
  print(paste0("Prior variance of s: ", s_eb/small_val  ))
  
  # Fit the model
  for (var_fct in vars_fct_NegBinBB_bayes){
    
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
  save(list_prior_fit_NegBinBB, vars_fct_NegBinBB_bayes, file =  "R_script_paper/fullybayes_Plants_fit_estimate_NegBinBBcompetitor.RData")
  
}


# We also consider classicBB, classicIBP and PoissonBB as competitor for BayesFactor
if (!file.exists("R_script_paper/fullybayes_Plants_fit_estimate_classics_and_PoissonBB.RData")) {
  
  # A) Mixtures of BBs (classicBB and PoissonBB)
  
  # MCMC setting for both 
  mcmcparams_both <- list(tau = 0.1, 
                          S = 5*10^4, n_burnin = 5*10^3, thin = 2)
  mcmcparams_obj_both <- mcmcparameters(model = "classicBB", mcmcparams = mcmcparams_both)
  
  # 1) classicBB
  # Initialization
  init_classicBB <- list(alpha_0 = - 1, s_0 = 15)
  init_obj_classicBB <- initialization(model = "classicBB", init = init_classicBB )
  
  # Prior: EB estimates
  small_val <- 10^(-4) # 10^(-3) - value in the first manuscript
  alpha_eb <- eb_EFPF_fit_PoissonBB$alpha
  theta_eb <- eb_EFPF_fit_PoissonBB$theta
  N_eb <- eb_EFPF_fit_PoissonBB$lambda
  
  s_eb <- alpha_eb + theta_eb
  
  print(paste0("Prior variance of -alpha: ", - alpha_eb/small_val  ))
  print(paste0("Prior variance of s: ", s_eb/small_val  ))
  
  # Prior: set hyperparameters
  hyper_classicBB <- list(a_alpha = - alpha_eb*small_val, b_alpha = small_val,
                          a_s = s_eb*small_val , b_s = small_val,
                          N = N_eb)
  prior_obj_classicBB <- prior(model = "classicBB", hyper = hyper_classicBB)
  
  # Fit the model
  prior_fit_classicBB <- GibbsFA(feature_matrix = data_mat,
                                 model = "classicBB", 
                                 prior = prior_obj_classicBB,
                                 initialization = init_obj_classicBB,
                                 mcmcparams = mcmcparams_obj_both)
  
  
  # 2) PoissonBB
  # Initialization
  init_PoissonBB <- list(alpha_0 = - 0.001, s_0 = 2)
  init_obj_PoissonBB <- initialization(model = "PoissonBB", init = init_PoissonBB )
  
  # Prior: EB estimates
  small_val <- 10^(-4) # 10^(-3) - value in the first manuscript
  alpha_eb <- eb_EFPF_fit_PoissonBB$alpha
  theta_eb <- eb_EFPF_fit_PoissonBB$theta
  lambda_eb <- eb_EFPF_fit_PoissonBB$lambda
  
  s_eb <- alpha_eb + theta_eb
  
  print(paste0("Prior variance of -alpha: ", - alpha_eb/small_val  ))
  print(paste0("Prior variance of s: ", s_eb/small_val  ))
  
  # Prior: set hyperparameters
  hyper_PoissonBB <- list(a_alpha = - alpha_eb*small_val, b_alpha = small_val,
                          a_s = s_eb*small_val , b_s = small_val,
                          lambda = lambda_eb)
  prior_obj_PoissonBB <- prior(model = "PoissonBB", hyper = hyper_PoissonBB)
  
  # Fit the model
  prior_fit_PoissonBB <- GibbsFA(feature_matrix = data_mat,
                                 model = "PoissonBB", 
                                 prior = prior_obj_PoissonBB,
                                 initialization = init_obj_PoissonBB,
                                 mcmcparams = mcmcparams_obj_both)
  
  
  # B) Mixtures of IBPs (classicIBP)
  
  # MCMC setting
  mcmcparams_classicIBP <- list(sigq_alpha = 0.1, sigq_s = 0.1, 
                                S = 5*10^4, n_burnin = 5*10^3, thin = 2)
  mcmcparams_obj_classicIBP <- mcmcparameters(model = "classicIBP", mcmcparams = mcmcparams_classicIBP)
  
  # Initialization
  init_classicIBP <- list(alpha_0 = 0.1, s_0 = 2)
  init_obj_classicIBP <- initialization(model = "classicIBP", init = init_classicIBP )
  
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
  
  # Prior: set hyperparameters
  gam_eb <- list_eb_EFPF_fit_GammaIBP[[1]]$gam
  
  # Hyperparameters elicitation 
  hyper_classicIBP <- list(gam = gam_eb,
                           a_alpha = small_val_alpha, b_alpha = t_eb*small_val_alpha,
                           a_s = s_eb*small_val_s , b_s = small_val_s)
  prior_obj_classicIBP <- prior(model = "classicIBP", hyper = hyper_classicIBP)
  
  
  prior_fit_classicIBP <- 
    GibbsFA(feature_matrix = data_mat,
            model = "classicIBP", 
            prior = prior_obj_classicIBP,
            initialization = init_obj_classicIBP,
            mcmcparams = mcmcparams_obj_classicIBP)
  
  
  
  # Save the entire workspace related to the type just performed
  save(prior_fit_classicBB, 
       prior_fit_PoissonBB, 
       prior_fit_classicIBP, 
       file =  "R_script_paper/fullybayes_Plants_fit_estimate_classics_and_PoissonBB.RData")
  
}




# Load the Work space
load("R_script_paper/fullybayes_Plants_fit_estimate_GammaIBP.RData")
load("R_script_paper/fullybayes_Plants_fit_estimate_classics_and_PoissonBB.RData")
load("R_script_paper/fullybayes_Plants_fit_estimate_NegBinBBcompetitor.RData")
if (!all(vars_GammaIBP == vars_GammaIBP_bayes)){
  stop("EB and FullyBayes use different prior variances for GammaIBP models")
}
if (!all(vars_fct_NegBinBB == vars_fct_NegBinBB_bayes)){
  stop("EB and FullyBayes use different prior variances for NegBinBB models")
}


## Convergence checks --------
library(ggmcmc)
library(coda)

# GammaIBP + prior
params_prior_GammaIBP <- list_prior_fit_GammaIBP[[paste0("var.", vars_GammaIBP[2])]][c("a_chain","b_chain", "alpha_chain", "theta_chain")]
params_prior_GammaIBP_df <- as.data.frame(do.call(cbind, params_prior_GammaIBP))

samples_GammaIBP <- mcmc.list(mcmc(params_prior_GammaIBP_df))
samples_ggs_GammaIBP <- ggs(samples_GammaIBP, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_GammaIBP)

effectiveSize(params_prior_GammaIBP_df)

# classicIBP + prior
params_prior_classicIBP <- prior_fit_classicIBP[c("alpha_chain", "theta_chain")]
params_prior_classicIBP_df <- as.data.frame(do.call(cbind, params_prior_classicIBP))

samples_classicIBP <- mcmc.list(mcmc(params_prior_classicIBP_df))
samples_ggs_classicIBP <- ggs(samples_classicIBP, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_classicIBP)

effectiveSize(params_prior_classicIBP_df)

# NegBinBB + prior
params_prior_NegBinBB <- list_prior_fit_NegBinBB[[paste0("var_fct.", vars_fct_NegBinBB[1])]][c("n0_chain","mu0_chain", "alpha_chain", "theta_chain")]
params_prior_NegBinBB_df <- as.data.frame(do.call(cbind, params_prior_NegBinBB))

samples_NegBinBB <- mcmc.list(mcmc(params_prior_NegBinBB_df))
samples_ggs_NegBinBB <- ggs(samples_NegBinBB, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_NegBinBB)

effectiveSize(params_prior_NegBinBB_df)

# classicBB + prior
params_prior_classicBB <- prior_fit_classicBB[c("alpha_chain", "theta_chain")]
params_prior_classicBB_df <- as.data.frame(do.call(cbind, params_prior_classicBB))

samples_classicBB <- mcmc.list(mcmc(params_prior_classicBB_df))
samples_ggs_classicBB <- ggs(samples_classicBB, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_classicBB)

effectiveSize(params_prior_classicBB_df)

# PoissonBB + prior
params_prior_PoissonBB <- prior_fit_PoissonBB[c("alpha_chain", "theta_chain")]
params_prior_PoissonBB_df <- as.data.frame(do.call(cbind, params_prior_PoissonBB))

samples_PoissonBB <- mcmc.list(mcmc(params_prior_PoissonBB_df))
samples_ggs_PoissonBB <- ggs(samples_PoissonBB, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_PoissonBB)

effectiveSize(params_prior_PoissonBB_df)



## Formal model-checking via BAYES FACTOR -------
log_marginal_like_object_list <- vector("list", length = 0)

for (var_fct_NegBinBB in vars_fct_NegBinBB){
  log_marginal_like_object_list[[paste0("NegBinBB.var_fct.", var_fct_NegBinBB)]] <- compute_log_marginal_likelihood_bridge(
    list_prior_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]]
  )
}
log_marginal_like_object_list[["classicBB"]] <- compute_log_marginal_likelihood_bridge(
  prior_fit_classicBB)
log_marginal_like_object_list[["PoissonBB"]] <- compute_log_marginal_likelihood_bridge(
  prior_fit_PoissonBB)
for (var_GammaIBP in vars_GammaIBP){
  log_marginal_like_object_list[[paste0("GammaIBP.var.", var_GammaIBP)]] <- compute_log_marginal_likelihood_bridge(
    list_prior_fit_GammaIBP[[paste0("var.", var_GammaIBP)]]
  )
}
log_marginal_like_object_list[["classicIBP"]] <- compute_log_marginal_likelihood_bridge(
  prior_fit_classicIBP)


print(log_marginal_like_object_list) 


log_marginal_like_list <- lapply(log_marginal_like_object_list, function(x)
  x$logml)
log_marginal_like_df <- data.frame(
  name = names(log_marginal_like_list),
  value = unlist(log_marginal_like_list),
  row.names = NULL
)
print(log_marginal_like_df)
#write.csv(log_marginal_like_df, file = "log_marginal_like_df.csv", row.names = FALSE)


### Tables of log-Bayes Factors
# Get all unique unordered combinations of names
name_combos <- combn(names(log_marginal_like_list), 2)

# Compute logBF and store results
logBF_names <- apply(name_combos, 2, function(x) paste(x[1], "vs", x[2], sep = "_"))
logBF_values <- apply(name_combos, 2, function(x) 
  log_marginal_like_list[[x[1]]] - log_marginal_like_list[[x[2]]])

# Assemble into data frame
logBF_df <- data.frame(Models = logBF_names, logBF = logBF_values)
print(logBF_df)



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



# HELD-OUT analysis with EFPF approach ------
frac_train <- 0.5
n_train <- ceiling(n*frac_train)
n_test <- n - n_train

set.seed(1234)
ind_train <- sample(1:n, n_train, replace = F)
data_mat_train <- data_mat[ind_train, ] 
data_mat_test <- data_mat[-ind_train, ]
data_mat_full <- rbind(data_mat_train, data_mat_test)

data_list_train <- convert_features_list(data_mat_train)  
feature_labels_train <- unique(unlist(data_list_train))
Kn_train <- length(feature_labels_train)

data_list_test <- convert_features_list(data_mat_test)  


## Train the models on training data -----

# Choices of variances
vars_fct_NegBinBB <- c(2,10) 
vars_GammaIBP <- c(10, 100) 

# Initial parameters for optimization
eb_init_BB <- list(alpha = -10, s = 100, Nhat_prime = 100)
eb_known_BB <- list()

eb_init_IBP <- list(alpha = 0.5, s = 10, Gamma = 10)
eb_known_IBP <- list()

eb_params_obj_BB <- eb_params(model = "BB", 
                              init = eb_init_BB, known = eb_known_BB )
eb_params_obj_IBP <- eb_params(model = "IBP", 
                               init = eb_init_IBP, known = eb_known_IBP )

# PoissonBB
eb_EFPF_fit_PoissonBB <- GibbsFA_eb(feature_matrix = data_mat_train, 
                                    model = "PoissonBB_eb", 
                                    type = "EFPF",
                                    eb_params =  eb_params_obj_BB)

# NegBinBB
list_eb_EFPF_fit_NegBinBB <-  vector(mode = "list", length = length(vars_fct_NegBinBB))
names(list_eb_EFPF_fit_NegBinBB) <- paste0("var_fct.", vars_fct_NegBinBB)

for (var_fct_NegBinBB in vars_fct_NegBinBB){
  
  list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]] <- 
    GibbsFA_eb(feature_matrix = data_mat_train,
               model = "NegBinBB_eb", type = "EFPF",
               eb_params =  eb_params_obj_BB, 
               var_fct = var_fct_NegBinBB)
  
}

# GammaIBP
list_eb_EFPF_fit_GammaIBP <-  vector(mode = "list", length = length(vars_GammaIBP))
names(list_eb_EFPF_fit_GammaIBP) <- paste0("var.", vars_GammaIBP)

for (var_GammaIBP in vars_GammaIBP){
  
  list_eb_EFPF_fit_GammaIBP[[paste0("var.", var_GammaIBP)]] <-
    GibbsFA_eb(feature_matrix = data_mat_train,
               model = "GammaIBP_eb", type = "EFPF",
               eb_params =  eb_params_obj_IBP,
               var_GammaIBP = var_GammaIBP)
  
}



## OPTION 1) Evaluate trained models on test data (and some reshuffles) -------

### 1) Compute prediction on number of new features: this just depends on the training data -------

# PoissonBB
n_new_EFPF_PoissonBB <- tibble(lambda = unname(unlist( 
  extrapolation(object = eb_EFPF_fit_PoissonBB, M = n_test, seed = seed)$lambda_post))) %>%
  mutate(lb = qpois(0.025, lambda, lower.tail = TRUE, log.p = FALSE),
         ub = qpois(0.975, lambda, lower.tail = TRUE, log.p = FALSE)) %>%
  rename(mean = lambda) %>%
  add_column(m = 1:n_test, Model = "Poisson BB") %>%
  select(mean, lb, ub, m, Model)


# NegBinBB
list_n_new_EFPF_NegBinBB <- vector(mode = "list", length = length(vars_fct_NegBinBB))
names(list_n_new_EFPF_NegBinBB) <- paste0("var_fct.", vars_fct_NegBinBB)

for (var_fct_NegBinBB in vars_fct_NegBinBB){
  
  eb_EFPF_NegBinBB_var <- list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]]
  
  list_n_new_EFPF_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]]  <- tibble(mu0 = unname(unlist(
    extrapolation(object = eb_EFPF_NegBinBB_var, M = n_test, seed = seed)$mu0_post )),
    n0 = unname(unlist( extrapolation(object = eb_EFPF_NegBinBB_var, M = n_test, seed = seed)$n0_post ))) %>%
    mutate(p = 1/(mu0/n0 + 1),
           lb = qnbinom(0.025, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE),
           ub = qnbinom(0.975, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE)) %>%
    rename(mean = mu0) %>%
    add_column(m = 1:n_test) %>%
    select(mean, lb, ub, m)
  
}


# GammaIBP
list_n_new_EFPF_GammaIBP <- vector(mode = "list", length = length(vars_GammaIBP))
names(list_n_new_EFPF_GammaIBP) <- paste0("var.", vars_GammaIBP)

for (var_GammaIBP in vars_GammaIBP){
  
  eb_EFPF_GammaIBP_var <- list_eb_EFPF_fit_GammaIBP[[paste0("var.", var_GammaIBP)]]
  
  list_n_new_EFPF_GammaIBP[[paste0("var.", var_GammaIBP)]] <- tibble(mu0 = unname(unlist(
    extrapolation(object = eb_EFPF_GammaIBP_var, M = n_test, seed = seed)$mu0_post )),
    n0 = unname(unlist( extrapolation(object = eb_EFPF_GammaIBP_var, M = n_test, seed = seed)$n0_post ))) %>%
    mutate(p = 1/(mu0/n0 + 1),
           lb = qnbinom(0.025, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE),
           ub = qnbinom(0.975, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE)) %>%
    rename(mean = mu0) %>%
    add_column(m = 1:n_test) %>%
    select(mean, lb, ub, m)
  
}


### 2) Compute observed number of new features from 1 to n_test, for each reshuffle -----
n_reshuffles <- 10

# Data-frame to store the number of hitherto unseen features in the test 
df_n_new_test <- matrix(nrow = n_reshuffles, ncol = n_test)

# Loop over the different reshuffles
for (d in 1:n_reshuffles){
  
  set.seed(123 + d)
  data_list_test_reshuffled <- sample(data_list_test)
  
  for (m in 1:n_test){
    data_list_test_reshuffled_first_m <- data_list_test_reshuffled[1:m]
    feature_labels_test_first_m <- unique(unlist(data_list_test_reshuffled_first_m))
    
    df_n_new_test[d, m] <- length(setdiff(feature_labels_test_first_m, feature_labels_train))
    
  }
  
}

# Summarize observed information in confidence intervals of reshuffles, for each m
# Function to compute empirical CI
get_empirical_ci <- function(x, probs = c(0.025, 0.975)) {
  mean_x <- mean(x)
  quantiles <- quantile(x, probs = probs)
  c(mean = mean_x, lb = unname(quantiles[1]), ub = unname(quantiles[2]))
}

n_new_test_mean_ci <- as_tibble(t(apply(df_n_new_test, 2, get_empirical_ci))) %>%
  mutate(m = row_number())



## OPTION 2) Plot extrapolation as in simulations  -----

### 1) Compute extrapolation: this just depends on the training data -------

# PoissonBB
extr_EFPF_PoissonBB_df <- tibble(lambda = unname(unlist( 
  extrapolation(object = eb_EFPF_fit_PoissonBB, M = n_test, seed = seed)$lambda_post))) %>%
  mutate(lb = qpois(0.025, lambda, lower.tail = TRUE, log.p = FALSE),
         ub = qpois(0.975, lambda, lower.tail = TRUE, log.p = FALSE)) %>%
  rename(mean = lambda) %>%
  add_row(mean = 0, lb = 0, ub = 0) %>%
  mutate(mean = mean + Kn_train, lb = lb + Kn_train, ub = ub + Kn_train) %>%
  add_column(x = c((n_train+1):(n_train + n_test), n_train), Model = "Poisson BB") %>%
  select(mean, lb, ub, x, Model)

extr_EFPF_PoissonBB_df$x <- as.integer(extr_EFPF_PoissonBB_df$x)


# NegBinBB
extr_EFPF_NegBinBB_df <- tibble(mean = numeric(), 
                                lb = numeric(), ub = numeric(),
                                x = integer(), Model = character())

for (var_fct_NegBinBB in vars_fct_NegBinBB){
  
  eb_EFPF_NegBinBB_var <- list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]]
  
  extr_EFPF_NegBinBB_df_var <- tibble(mu0 = unname(unlist(
    extrapolation(object = eb_EFPF_NegBinBB_var, M = n_test, seed = seed)$mu0_post )),
    n0 = unname(unlist( extrapolation(object = eb_EFPF_NegBinBB_var, M = n_test, seed = seed)$n0_post ))) %>%
    mutate(p = 1/(mu0/n0 + 1),
           lb = qnbinom(0.025, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE),
           ub = qnbinom(0.975, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE)) %>%
    rename(mean = mu0) %>%
    add_row(mean = 0, lb = 0, ub = 0) %>%
    mutate(mean = mean + Kn_train, lb = lb + Kn_train, ub = ub + Kn_train) %>%
    add_column(x = c((n_train+1):(n_train + n_test), n_train), Model = paste0("NegBinomial BB x", var_fct_NegBinBB)) %>%
    select(mean, lb, ub, x, Model)
  
  extr_EFPF_NegBinBB_df_var$x <- as.integer(extr_EFPF_NegBinBB_df_var$x)
  
  extr_EFPF_NegBinBB_df <- bind_rows(extr_EFPF_NegBinBB_df, 
                                     extr_EFPF_NegBinBB_df_var)
  
  
}


# GammaIBP
extr_EFPF_GammaIBP_df <- tibble(mean = numeric(), 
                                lb = numeric(), ub = numeric(),
                                x = integer(), Model = character())

for (var_GammaIBP in vars_GammaIBP){
  
  eb_EFPF_GammaIBP_var <- list_eb_EFPF_fit_GammaIBP[[paste0("var.", var_GammaIBP)]]
  
  extr_EFPF_GammaIBP_df_var <- tibble(mu0 = unname(unlist(
    extrapolation(object = eb_EFPF_GammaIBP_var, M = n_test, seed = seed)$mu0_post )),
    n0 = unname(unlist( extrapolation(object = eb_EFPF_GammaIBP_var, M = n_test, seed = seed)$n0_post ))) %>%
    mutate(p = 1/(mu0/n0 + 1),
           lb = qnbinom(0.025, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE),
           ub = qnbinom(0.975, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE)) %>%
    rename(mean = mu0) %>%
    add_row(mean = 0, lb = 0, ub = 0) %>%
    mutate(mean = mean + Kn_train, lb = lb + Kn_train, ub = ub + Kn_train) %>%
    add_column(x = c((n_train+1):(n_train + n_test), n_train), Model = paste0("GammaIBP, var = ", var_GammaIBP)) %>%
    select(mean, lb, ub, x, Model)
  
  
  extr_EFPF_GammaIBP_df_var$x <- as.integer(extr_EFPF_GammaIBP_df_var$x)
  
  extr_EFPF_GammaIBP_df <- bind_rows(extr_EFPF_GammaIBP_df, 
                                     extr_EFPF_GammaIBP_df_var)
  
}



### 2) Accumulation curve for training and test set -----

# Extract accumulation curve of the observed sample (or average accumulation)
accum_df <- tibble(n_feat = c(0,rarefaction(data_mat_full, n_reorderings = 1)),
                   type = c(rep("train", n_train +1), rep("test", n_test)),
                   x = 0:(n_train + n_test))

accum_df_train <- accum_df %>%
  filter(type == "train")
accum_df_test <- accum_df %>%
  filter(type == "test")


### 3) Plot -----

extr_all_df <- bind_rows(extr_EFPF_PoissonBB_df, 
                         extr_EFPF_NegBinBB_df,
                         extr_EFPF_GammaIBP_df) #%>% filter(Model == "NegBinomial BB x10")
#%>% filter(Model %in% c("Poisson BB", paste0("NegBinomial BB x", vars_fct_NegBinBB)))

extr_all_df$Model <- factor(extr_all_df$Model, 
                            levels = c("Poisson BB", 
                                       paste0("NegBinomial BB x", vars_fct_NegBinBB),
                                       paste0("GammaIBP, var = ", vars_GammaIBP)))

get_blue_tones <- function(n) {
  all_blues <- c(
    "blue4", "dodgerblue4", "dodgerblue2", "lightskyblue", "lightblue")
  
  if (n > length(all_blues)) {
    stop("Requested number exceeds available blue tones.")
  }
  
  return(all_blues[1:n])
}

get_red_tones <- function(n) {
  all_reds <- c(
    "darkred", "tomato3", "lightsalmon")
  
  if (n > length(all_reds)) {
    stop("Requested number exceeds available red tones.")
  }
  
  return(all_reds[1:n])
}

models_name <- unique(extr_all_df$Model)
number_bbs_models <- sum(grepl("^Pois", models_name)) + sum(grepl("^Neg", models_name))
number_ibps_models <- sum(grepl("^Gam", models_name))
colors_name <- c(get_blue_tones(number_bbs_models), get_red_tones(number_ibps_models)) 
color_dict <- setNames(colors_name, models_name)

plot_ribbons <- ggplot() +
  geom_point(data = accum_df_train, aes(x = x, y = n_feat),
             color="black", shape = 19, size = 0.5) +
  geom_point( data = accum_df_test, aes(x = x, y = n_feat),
              color="black", shape = 19, size = 0.5) +
  geom_ribbon(data = extr_all_df, aes(x = x, ymin = lb, ymax = ub, fill = Model), color = NA, alpha = 0.4) +
  scale_fill_manual(values = color_dict) +
  # facet_grid(~ n_train_latex,
  #            labeller = label_parsed,
  #            scales = "free_x")  +
  #geom_vline(data = temp, mapping =  aes(xintercept = xvalues) , linetype = "dashed", color = "grey") +
  geom_vline(xintercept = n_train, linetype = "dashed", color = "grey") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) +
  scale_color_tableau()

plot_ribbons
