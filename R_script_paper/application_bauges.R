rm(list=ls())

library(ProductFormFA)
library(tidyverse)
library(scales)
library(ggpubr)
library(ggthemes)

source("R_script_paper/Routine_Chao.R")
source("R_script_paper/utils.R")

# Fit models -----

data <- read.csv(file = "R_script_paper/Bauges_data.csv", header = TRUE,
                 row.names="X")

data <- as.matrix(data)
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
plot(1:n, rarefaction(data_mat, n_reorderings = 10) )


# 2) Run the models on the data 

var_fct_NegBinBB <- 20
var_GammaIBP <- 10

eb_init_BB <- list(alpha = -10, s = 100, Nhat_prime = 200)
eb_known_BB <- list()

eb_init_IBP <- list(alpha = 0.5, s = 1, Gamma = 10)
eb_known_IBP <- list()

eb_params_obj_BB <- eb_params(model = "BB", 
                              init = eb_init_BB, known = eb_known_BB )
eb_params_obj_IBP <- eb_params(model = "IBP", 
                               init = eb_init_IBP, known = eb_known_IBP )

# PoissonBB
eb_EFPF_PoissonBB_fit <- GibbsFA_eb(feature_matrix = data_mat, 
                                    model = "PoissonBB", 
                                    type = "EFPF",
                                    eb_params =  eb_params_obj_BB)

# NegBinBB
eb_EFPF_NegBinBB_fit <- GibbsFA_eb(feature_matrix = data_mat,
                                   model = "NegBinBB", type = "EFPF",
                                   eb_params =  eb_params_obj_BB, 
                                   var_fct = var_fct_NegBinBB)


# GammaIBP
eb_EFPF_GammaIBP_fit <- GibbsFA_eb(feature_matrix = data_mat,
                                   model = "GammaIBP", type = "EFPF",
                                   eb_params =  eb_params_obj_IBP,
                                   var_GammaIBP = var_GammaIBP)




# Model Checking ----

# 0.B) Check on rarefaction
n_rare <- n
eb_EFPF_fit_PoissonBB_rare <- eb_EFPF_PoissonBB_fit
eb_EFPF_fit_NegBinBB_rare <- eb_EFPF_NegBinBB_fit
eb_EFPF_fit_GammaIBP_rare <- eb_EFPF_GammaIBP_fit

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

df_rare <- rbind(rare_EFPF_PoissonBB,
                 rare_EFPF_NegBinBB,
                 rare_EFPF_GammaIBP)

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
ggsave(filename = "R_script_paper/Paper_plots/rarefaction_Bauges_eb_EFPF.pdf", width = 8, height = 3.5, dpi = 300, units = "in", device='pdf')


# 0.c) Check on K_n_r
n_knr <- n
eb_EFPF_fit_PoissonBB_knr <- eb_EFPF_PoissonBB_fit
eb_EFPF_fit_NegBinBB_knr <- eb_EFPF_NegBinBB_fit
eb_EFPF_fit_GammaIBP_knr <- eb_EFPF_GammaIBP_fit

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

df_K_n_r <- rbind(K_n_r_EFPF_PoissonBB,
                  K_n_r_EFPF_NegBinBB,
                  K_n_r_EFPF_GammaIBP)

df_K_n_r$Model <- factor(df_K_n_r$Model)

r_positive <- observed_K_n_r %>%
  filter(k_n_r > 0) %>%
  select(r) %>%
  filter(r < 100)

df_K_n_r_plot <- df_K_n_r %>%
  filter(r %in% c(r_positive$r))

observed_K_n_r_plot <- observed_K_n_r %>%
  filter(r %in% c(r_positive$r))

ggplot(df_K_n_r_plot, aes(x = r, y = means, color = Model)) +
  geom_line( linetype = "dashed") +
  geom_point( data = observed_K_n_r_plot, aes(x = r, y = k_n_r), color="black", size = 1) +
  scale_y_log10() +
  xlab("r") + ylab(expression(m[r])) + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1)
ggsave(filename = paste0("R_script_paper/Paper_plots/knr_", type, "_eb_EFPF.pdf"), width = 5, height = 4.5, dpi = 300, units = "in", device='pdf')



# Prediction ----

## Richness estimation -----

params_richness_EFPF_PoissonBB <- tibble( lambda_prime = 
                                            total_richness(eb_EFPF_PoissonBB_fit)$lambda_post) %>%
  add_column(Model = "PoissonBB") %>%
  mutate(lb = qpois(0.025, lambda_prime, lower.tail = TRUE, log.p = FALSE),
         ub = qpois(0.975, lambda_prime, lower.tail = TRUE, log.p = FALSE) )


params_richness_EFPF_NegBinBB <- tibble( n0_prime = 
                                           total_richness(eb_EFPF_NegBinBB_fit)$n0_post,
                                         mu0_prime = 
                                           total_richness(eb_EFPF_NegBinBB_fit)$mu0_post) %>%
  add_column(Model = "NegBinBB") %>%
  mutate(p_prime = 1/(mu0_prime/n0_prime + 1),
         lb = qnbinom(0.025, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE),
         ub = qnbinom(0.975, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE) )



bounds <- tibble( lb = min(params_richness_EFPF_PoissonBB$lb + Kn, 
                           params_richness_EFPF_NegBinBB$lb + Kn),
                  ub = max(params_richness_EFPF_PoissonBB$ub + Kn,
                           params_richness_EFPF_NegBinBB$ub + Kn))


dens_richness_PoissonBB <- tibble( x = bounds$lb: bounds$ub) %>%
  mutate( y = dpois(x - Kn , lambda = params_richness_EFPF_PoissonBB$lambda_prime)) %>%
  add_column(Model = "PoissonBB")

dens_richness_NegBinBB <- tibble( x = bounds$lb : bounds$ub) %>%
  mutate( y = dnbinom(x - Kn, size = params_richness_EFPF_NegBinBB$n0_prime, 
                      prob = params_richness_EFPF_NegBinBB$p_prime)) %>%
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
ggsave(filename = paste0("R_script_paper/Paper_plots/richness_", type, "_eb_EFPF.pdf"), width = 3.8, height = 3.8, dpi = 300, units = "in", device='pdf')


## Extrapolation (EFPF version) -----

# Extract accumulation curve of the observed sample (or average accumulation)
M = 1000

accum_df <- tibble( x = 0:n,
                    n_feat = c(0,rarefaction(data_mat, n_reorderings = 20)))


extr_PoissonBB_df <- tibble(lambda = unname(unlist( 
  extrapolation(object = eb_EFPF_PoissonBB_fit, M = M, seed = seed)$lambda_post)),
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
  extrapolation(object = eb_EFPF_NegBinBB_fit, M = M, seed = seed)$mu0_post )),
  n0 = unname(unlist( extrapolation(object = eb_EFPF_NegBinBB_fit, M = M, seed = seed)$n0_post )),
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
  extrapolation(object = eb_EFPF_GammaIBP_fit, M = M, seed = seed)$mu0_post )),
  n0 = unname(unlist( extrapolation(object = eb_EFPF_GammaIBP_fit, M = M, seed = seed)$n0_post )),
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
  geom_ribbon(aes(ymin = lb, ymax = ub, color = Model), linewidth = 0.8, alpha = 0.1) +
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
ggsave(filename = paste0("R_script_paper/Paper_plots/extr_", type, "_eb_EFPF.pdf"), width = 3.8, height = 3.8, dpi = 300, units = "in", device='pdf')







# Prior approach - to be done ----------

# Empirical estimate of E(N) is obtained by Chiu
Nbar <- beta_binomial_estimator(data_mat)

# 2.1) Beta-Bernoulli with Poisson(lambda) mixture

# Initialization and MCMC setting 
init_PoissonBB <- list(alpha_0 = -1, s_0 = 1)
init_obj_PoissonBB <- initialization(model = "PoissonBB", init = init_PoissonBB )
mcmcparams_PoissonBB <- list(tau = 0.1, S = 3*10^4, n_burnin = 5*10^3, thin = 2)
mcmcparams_obj_PoissonBB <- mcmcparameters(model = "PoissonBB", mcmcparams = mcmcparams_PoissonBB)

# Hyperparameters elicitation 
hyper_PoissonBB <- list(a_alpha = 1, b_alpha = 0.1,
                        a_s = 2, b_s = 0.2,
                        lambda = Nbar)
prior_obj_PoissonBB <- prior(model = "PoissonBB", hyper = hyper_PoissonBB) 

# Fit the model
PoissonBB_fit <- GibbsFA(feature_matrix = data_mat, 
                         model = "PoissonBB", 
                         prior = prior_obj_PoissonBB, 
                         initialization = init_obj_PoissonBB, 
                         mcmcparams = mcmcparams_obj_PoissonBB)

# 2.1.B) Beta-Bernoulli with Poisson(lambda) mixture - EB version

# EB parameters 
eb_init <- list(alpha = -1, s = 1)#, lambda = 100)
eb_known <- list(lambda = Nbar)
eb_params_obj <- eb_params(model = "PoissonBB", init = eb_init, known = eb_known )

# Fit the model with EB
eb_PoissonBB_fit <- GibbsFA_eb(feature_matrix = data_mat, 
                               model = "PoissonBB", 
                               eb_params =  eb_params_obj)




# 2.2) Beta-Bernoulli with NB(n0, mu0) mixture

# Initialization and MCMC setting
init_NegBinBB <- list(alpha_0 = -1, s_0 = 1)
init_obj_NegBinBB <- initialization(model = "NegBinBB", init = init_NegBinBB )
mcmcparams_NegBinBB <- list(tau = 0.1, S = 3*10^4, n_burnin = 5*10^3, thin = 2)
mcmcparams_obj_NegBinBB <- mcmcparameters(model = "NegBinBB", mcmcparams = mcmcparams_NegBinBB)

# Hyperparameters elicitation 
c_fr <- 10

hyper_NegBinBB <- list(a_alpha = 1, b_alpha = 0.1,
                       a_s = 2, b_s = 0.2,
                       n0 = Nbar/(c_fr - 1), # n0, mu0 are set s.t. E(N) = Nbar, Var(N) = c_fr*E(N)
                       mu0 = 1/c_fr)
prior_obj_NegBinBB <- prior(model = "NegBinBB", hyper = hyper_NegBinBB) 


# Fit the model
NegBinBB_fit <- GibbsFA(feature_matrix = data_mat, 
                        model = "NegBinBB", 
                        prior = prior_obj_NegBinBB, 
                        initialization = init_obj_NegBinBB, 
                        mcmcparams = mcmcparams_obj_NegBinBB)

# 2.2.B) Beta-Bernoulli with NB(n0,mu0) mixture - EB version

# EB parameters 
eb_init <- list(alpha = -1, s = 1)
c <- 10
eb_known <- list( mu0 = Nbar, n0 = Nbar/(c-1))
eb_params_obj <- eb_params(model = "NegBinBB", init = eb_init, known = eb_known )

# Fit the model with EB
eb_NegBinBB_fit <- GibbsFA_eb(feature_matrix = data_mat, 
                              model = "NegBinBB", 
                              eb_params =  eb_params_obj)


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








