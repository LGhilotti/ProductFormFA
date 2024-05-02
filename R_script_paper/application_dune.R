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

data <- as.matrix(dune > 0)

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


# 2.1) Beta-Bernoulli with Poisson(lambda) mixture - EB version

# EB parameters 
eb_init_PoissonBB <- list(alpha = -1000, s = 100, lambda = 100)
eb_known_PoissonBB <- list()
eb_params_obj_PoissonBB <- eb_params(model = "PoissonBB", 
                                     init = eb_init_PoissonBB, known = eb_known_PoissonBB )

eb_PoissonBB_fit <- GibbsFA_eb(feature_matrix = data_mat, 
                               model = "PoissonBB", 
                               eb_params =  eb_params_obj_PoissonBB)


eb_PoissonBB_fit

# 2.2) NegBinBB

# Initialization and known parameters
c_fr <- 5
eb_init_NegBinBB <- list(alpha = -100, s = 100, mu0 = 30)
eb_known_NegBinBB <- list(var_fct = c_fr)
eb_params_obj_NegBinBB <- eb_params(model = "NegBinBB",
                                    init = eb_init_NegBinBB,
                                    known = eb_known_NegBinBB)

eb_NegBinBB_fit <- GibbsFA_eb(feature_matrix = data_mat,
                              model = "NegBinBB",
                              eb_params =  eb_params_obj_NegBinBB)

eb_NegBinBB_fit


# 2.3) GammaIBP

# Initialization and known parameters
eb_init_GammaIBP <- list(alpha = 0.5, s = 1, a = 1, b = 1)
eb_known_GammaIBP <- list()
eb_params_obj_GammaIBP <- eb_params(model = "GammaIBP",
                                    init = eb_init_GammaIBP,
                                    known = eb_known_GammaIBP)

eb_GammaIBP_fit <- GibbsFA_eb(feature_matrix = data_mat,
                              model = "GammaIBP",
                              eb_params =  eb_params_obj_GammaIBP)

eb_GammaIBP_fit


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


# 0.B) Check on rarefaction
accum_df <- tibble( x = 0:n,
                    n_feat = c(0,rarefaction(data_mat[1:n,], n_reorderings = 20)))

rare_PoissonBB <- tibble( lambda_post = unname(unlist(
  rarefaction(object = eb_PoissonBB_fit, seed = seed)$lambda_post ))) %>%
  rename(means = lambda_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "PoissonBB",
             x = c(1:n,0))

rare_NegBinBB <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_NegBinBB_fit, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "NegBinBB",
             x = c(1:n,0))

rare_GammaIBP <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_GammaIBP_fit, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "GammaIBP",
             x = c(1:n,0)) 

df_rare <- rbind(rare_PoissonBB, rare_NegBinBB, rare_GammaIBP)
df_rare$Model <- factor(df_rare$Model, levels = c("PoissonBB", "NegBinBB", "GammaIBP"))

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
ggsave(filename = "R_script_paper/Paper_plots/rarefaction_dune_eb.pdf", width = 8, height = 3.5, dpi = 300, units = "in", device='pdf')






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








