#
### Functions for K^n_m plots (use n=0 to get prior K_n)
#

# mixtures of Beta-Bernoulli 

ci_kmn_bb <- function(alpha_bb, theta_bb, m, n = 0, k=0, N, lev){
  
  p_nm <- 1 - exp( lgamma(theta_bb + alpha_bb + n + (1:m)) - lgamma(theta_bb + alpha_bb + n) -
                     lgamma(theta_bb + n + (1:m)) + lgamma(theta_bb + n) )
  
  
  means_binom <- (N-k)*p_nm
  
  ubs <- qbinom(lev + (1 - lev) / 2, N-k, p_nm)
  lbs <- qbinom((1 - lev) / 2, N-k, p_nm) 
  
  return(list("means" = means_binom,"ubs" = ubs,"lbs" = lbs))
  
}


ci_kmn_poiss_bb <- function(alpha_bb, theta_bb, m, n = 0, lambda, lev){
  
  p_nm <- 1 - exp( lgamma(theta_bb + alpha_bb + n + (1:m)) - lgamma(theta_bb + alpha_bb + n) -
                     lgamma(theta_bb + n + (1:m)) + lgamma(theta_bb + n) )
  
  t <- exp( lgamma(theta_bb + alpha_bb + n) - lgamma(theta_bb + alpha_bb) -
              lgamma(theta_bb + n ) + lgamma(theta_bb) )
  
  poiss_pars <- p_nm*lambda*t
  
 
  ubs <- qpois(lev + (1 - lev) / 2, poiss_pars)
  lbs <- qpois((1 - lev) / 2, poiss_pars) 
  
  return(list("means" = poiss_pars,"ubs" = ubs,"lbs" = lbs))
  
}


ci_kmn_negbin_bb <- function(alpha_bb, theta_bb, m, n = 0, k = 0, n_0, p, lev){
  
  p_nm <- 1 - exp( lgamma(theta_bb + alpha_bb + n + (1:m)) - lgamma(theta_bb + alpha_bb + n) -
                     lgamma(theta_bb + n + (1:m)) + lgamma(theta_bb + n) )
  
  t <- exp( lgamma(theta_bb + alpha_bb + n) - lgamma(theta_bb + alpha_bb) -
              lgamma(theta_bb + n ) + lgamma(theta_bb) )
  
  pbar <- 1 - (p_nm*(1-p)*t)/(1 - (1-p_nm)*(1-p)*t )
  
  nb_first_par <- n_0 + k
  
  means_nb <- nb_first_par*(1-pbar)/pbar
  
  ubs <- qnbinom(lev + (1 - lev) / 2, nb_first_par, pbar)
  lbs <- qnbinom((1 - lev) / 2, nb_first_par, pbar) 
  
  return(list("means" = means_nb,"ubs" = ubs,"lbs" = lbs))
  
}

# mixtures of IBP

ci_kmn_ibp <- function(alpha_ibp, theta_ibp, m, n = 0, gam, lev){
  
  sum_m <- cumsum(exp( lgamma(theta_ibp + alpha_ibp + n + (1:m) - 1) - lgamma(theta_ibp + alpha_ibp) -
                     lgamma(theta_ibp + n + (1:m)) + lgamma(theta_ibp + 1) ))
  
  poiss_pars <- gam*sum_m
  
  
  ubs <- qpois(lev + (1 - lev) / 2, poiss_pars)
  lbs <- qpois((1 - lev) / 2, poiss_pars) 
  
  return(list("means" = poiss_pars,"ubs" = ubs,"lbs" = lbs))
  
}


ci_kmn_gamma_ibp <- function(alpha_ibp, theta_ibp, m, n = 0, k=0, a, b, lev){
  
  sum_m <- cumsum(exp( lgamma(theta_ibp + alpha_ibp + n + (1:m) - 1) - lgamma(theta_ibp + alpha_ibp) -
                         lgamma(theta_ibp + n + (1:m)) + lgamma(theta_ibp + 1) ))
  
  if (n == 0){
    g_n = 0
  } else {
    g_n <- sum(exp( lgamma(theta_ibp + alpha_ibp + (1:n) - 1) - lgamma(theta_ibp + alpha_ibp) -
                      lgamma(theta_ibp + (1:n)) + lgamma(theta_ibp + 1) ))
  }
  
  
  nb_first_par <- k + a
  nb_second_par <- (g_n + b)/(g_n + b + sum_m)
  
  means_nb <- nb_first_par*(1 - nb_second_par)/nb_second_par
  
  ubs <- qnbinom(lev + (1 - lev) / 2, nb_first_par, nb_second_par)
  lbs <- qnbinom((1 - lev) / 2, nb_first_par, nb_second_par) 
  
  return(list("means" = means_nb,"ubs" = ubs,"lbs" = lbs))
  
}

#
#### Kn plots ------------
#

library(tidyverse)
library(scales)
library(ggthemes)

lev <- 0.95
H <- 500 # we set all the mixtures of BB to have this expected value
lambda <- H
c_fr <- 5
n0_nb <- H/(c_fr - 1)
p_nb <- 1/c_fr

m <- 200
alpha_BB <- -1.2
theta_BB <- 10
alpha_IBP <- 0.8
theta_IBP <- 1
gam <- 10
a <- 100
b <- 10

ci_BB <- as_tibble( ci_kmn_bb(alpha_BB, theta_BB, m, n =0, k= 0, H, lev) ) %>%
  add_column(Model = "BB", t = 1:m)

ci_poiss_BB <- as_tibble( ci_kmn_poiss_bb(alpha_BB, theta_BB, m, n =0, lambda, lev) ) %>%
  add_column(Model = "Poisson", t = 1:m)

ci_negbin_BB <- as_tibble( ci_kmn_negbin_bb(alpha_BB, theta_BB, m, n =0,k=0, n_0 = n0_nb, p = p_nb, lev)) %>%
  add_column(Model = "NegBin", t = 1:m)

ci_IBP <- as_tibble( ci_kmn_ibp(alpha_IBP, theta_IBP, m, n =0, gam = gam, lev)) %>%
  add_column(Model = "3IBP", t = 1:m)

ci_gamma_IBP <- as_tibble( ci_kmn_gamma_ibp(alpha_IBP, theta_IBP, m, n =0, k=0, a = a, b = b, lev)) %>%
  add_column(Model = "Gamma", t = 1:m)

joint_bb <- rbind( ci_BB, ci_poiss_BB, ci_negbin_BB )
joint_ibp <- rbind( ci_IBP, ci_gamma_IBP )

ggplot(joint_bb, aes(t,means, color = Model)) +
  geom_line(linetype = "dashed", linewidth = 0.8) +
  geom_ribbon(aes(ymin = lbs, ymax = ubs), linewidth = 0.8, alpha = 0.1) +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_tableau() +
  theme(aspect.ratio = 1) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))



ggsave(filename = "Plots_paper/prior_kn_bb.pdf", width = 4, height = 4, dpi = 300, units = "in", device='pdf')

ggplot(joint_ibp, aes(t,means, color = Model)) +
  geom_line(linetype = "dashed", linewidth = 0.8) +
  geom_ribbon(aes(ymin = lbs, ymax = ubs), linewidth = 0.8, alpha = 0.1) +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_tableau() +
  theme(aspect.ratio = 1) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))



ggsave(filename = "Plots_paper/prior_kn_ibp.pdf", width = 4, height = 4, dpi = 300, units = "in", device='pdf')


#
#### K^n_m plots ------------
#

library(tidyverse)
library(scales)
library(ggthemes)
library(ProductFormFA)

m <- 200
n <- 40

# Observed sample <-> cts on the full sample
set.seed(1234)
buff_poiss_bb <- buffet_poiss_BB(alpha = - 1, theta = 8, n = n, lambda = 400)
plot_trajectory(buff_poiss_bb$features)
data_mat <- create_features_matrix(buff_poiss_bb$features)

obs_sample <- sapply(2:n, function(i) ncol(data_mat[1:i,colSums(data_mat[1:i,]) > 0])   )
obs_sample <- data.frame(t = 0:n, 
                         obs = c(0, sum(data_mat[1,]) , obs_sample))

k <- ncol(data_mat)

# Models parameter
lev <- 0.95
H <- 500 # we set all the mixtures of BB to have this expected value
lambda <- H
c_fr <- 5
n0_nb <- H/(c_fr - 1)
p_nb <- 1/c_fr

alpha_BB <- -1.2
theta_BB <- 10
alpha_IBP <- 0.8
theta_IBP <- 1
gam <- 10
a <- 100
b <- 10

ci_kmn_BB <- as_tibble( ci_kmn_bb(alpha_BB, theta_BB, m, n =n, k= k, N = H, lev) ) %>%
  add_column(Model = "BB", t = (n+1):(n+m) )

ci_kmn_poiss_BB <- as_tibble( ci_kmn_poiss_bb(alpha_BB, theta_BB, m, n =n, lambda, lev) ) %>%
  add_column(Model = "Poisson", t = (n+1):(n+m) )

ci_kmn_negbin_BB <- as_tibble( ci_kmn_negbin_bb(alpha_BB, theta_BB, m, n =n,k=k, n_0 = n0_nb, p = p_nb, lev)) %>%
  add_column(Model = "NegBin", t = (n+1):(n+m) )

ci_kmn_IBP <- as_tibble( ci_kmn_ibp(alpha_IBP, theta_IBP, m, n =n, gam = gam, lev)) %>%
  add_column(Model = "3IBP", t = (n+1):(n+m) )

ci_kmn_gamma_IBP <- as_tibble( ci_kmn_gamma_ibp(alpha_IBP, theta_IBP, m, n =n, k=k, a = a, b = b, lev)) %>%
  add_column(Model = "Gamma", t = (n+1):(n+m) )



joint_kmn_bb <- rbind( ci_kmn_BB, ci_kmn_poiss_BB, ci_kmn_negbin_BB )

joint_kmn_bb <- joint_kmn_bb %>%
  mutate(means = means + obs_sample[n+1,]$obs,
         lbs = lbs + obs_sample[n+1,]$obs,
         ubs = ubs + obs_sample[n+1,]$obs)


joint_kmn_ibp <- rbind( ci_kmn_IBP, ci_kmn_gamma_IBP )

joint_kmn_ibp <- joint_kmn_ibp %>%
  mutate(means = means + obs_sample[n+1,]$obs,
         lbs = lbs + obs_sample[n+1,]$obs,
         ubs = ubs + obs_sample[n+1,]$obs)


ggplot(joint_kmn_bb, aes(t,means, color = Model)) +
  geom_line(linetype = "dashed", linewidth = 0.8) +
  geom_ribbon(aes(ymin = lbs, ymax = ubs), linewidth = 0.8, alpha = 0.1) +
  geom_line(data = obs_sample, aes(t, obs), color="black", linetype="solid", linewidth=0.5) +
  geom_vline(xintercept = n , linetype = "dashed", color = "grey") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_tableau() +
  theme(aspect.ratio = 1) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) 


ggsave(filename = "Plots_paper/pred_kmn_bb.pdf", width = 4, height = 4, dpi = 300, units = "in", device='pdf')



ggplot(joint_kmn_ibp, aes(t,means, color = Model)) +
  geom_line(linetype = "dashed", linewidth = 0.8) +
  geom_ribbon(aes(ymin = lbs, ymax = ubs), linewidth = 0.8, alpha = 0.1) +
  geom_line(data = obs_sample, aes(t, obs), color="black", linetype="solid", linewidth=0.5) +
  geom_vline(xintercept = n , linetype = "dashed", color = "grey") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_tableau() +
  theme(aspect.ratio = 1) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))


ggsave(filename = "Plots_paper/pred_kmn_ibp.pdf", width = 4, height = 4, dpi = 300, units = "in", device='pdf')
