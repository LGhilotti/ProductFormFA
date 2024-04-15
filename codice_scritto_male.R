feature_fraction <- function(n, theta, alpha){
  1 - exp(lgamma(theta + alpha + n ) - lgamma(theta + alpha) - lgamma(theta + n) + lgamma(theta))
}
feature_fraction <- Vectorize(feature_fraction, "n")

gn <- function(n, theta, alpha){
  idx <- 1:n
  sum(exp(lgamma(theta + alpha + idx - 1) - lgamma(theta+alpha) + lgamma(theta + 1) - lgamma(theta + idx)))
}
gn <- Vectorize(gn, "n")

rarefaction_pois <- function(n, lambda, theta, alpha){
  sites <- 1:n
  richness <- lambda * feature_fraction(sites, theta, alpha)
  list(sites = sites, richness = richness)
}

richness_pois <- function(k, n, lambda, theta, alpha){
  richness <- k + lambda * (1 - feature_fraction(n, theta, alpha))
  richness
}

rarefaction_IBP <- function(n, theta, alpha, gamma){
  sites <- 1:n
  richness <- gamma * gn(sites, theta, alpha)
  list(sites = sites, richness = richness)
}

EB_EFPF_poiss_BB <- function(n, counts, pars_0) {
  # initialize s = theta + alpha
  pars_0[2] <- pars_0[2] + pars_0[1]
  
  # set constraints stricter so that function is limited
  res <- optim(
    par = pars_0, fn = neg_log_EFPF_poiss_BB_rep, n = n, counts = counts,
    method = "L-BFGS-B", lower = c(-Inf, 1e-5, 1e-5), upper = c(-1e-5, Inf, Inf)
  )
  
  sol <- res$par
  # convert to the theta parameter
  sol[2] <- sol[2] - sol[1]
  
  return(sol)
}


# Dataset
library(tidyverse)
fungi2016_mat <- read.csv("clusters.csv", sep = ";") %>% column_to_rownames(var = "CLUSTER")
fungi2016_mat[is.na(fungi2016_mat)] <- 0
#fungi2016_mat <- read_excel('mazziotta2016_application/mazz2016_data.xls', sheet = "Plants") %>% column_to_rownames(var = "Species")
fungi2016_mat <- t(fungi2016_mat)

# set.seed(012345)
# fungi2016_mat <- cbind(matrix(rbinom(50 * 10, size = 1, prob = 0.9), 50, 20),
#                        matrix(rbinom(50 * 10, size = 1, prob = 0.8), 50, 20),
#                        matrix(rbinom(50 * 10, size = 1, prob = 0.7), 50, 20),
#                        matrix(rbinom(50 * 10, size = 1, prob = 0.6), 50, 20),
#                        matrix(rbinom(50 * 10, size = 1, prob = 0.5), 50, 20),
#                        matrix(rbinom(50 * 10, size = 1, prob = 0.4), 50, 20),
#                        matrix(rbinom(50 * 10, size = 1, prob = 0.3), 50, 20),
#                        matrix(rbinom(50 * 10, size = 1, prob = 0.2), 50, 20),
#                        matrix(rbinom(50 * 100, size = 1, prob = 0.01), 50, 100))
fungi2016_mat <- fungi2016_mat[, colSums(fungi2016_mat) > 0]
n <- nrow(fungi2016_mat)
k <- ncol(fungi2016_mat)

hist(apply(fungi2016_mat, 2, mean))

set.seed(012345)
# This randomly choose an ordering
fungi2016_mat <-fungi2016_mat[sample.int(n, size = n, replace = F),]
# Conversion of the matrix into a list
fungi2016_list <- create_features_list(fungi2016_mat)

# Frequencies m_1,...,m_k
freq <- colSums(fungi2016_mat)
table(freq)

# Poisson case (finite N)
fit_pois <- EB_EFPF_poiss_BB(n = n, freq, pars_0 = c(-0.6, 1, 30))
alpha_hat <- fit_pois[1]
theta_hat <- fit_pois[2]
lambda_hat <- fit_pois[3]

alpha_hat
-alpha_hat / theta_hat
theta_hat

richness_pois(k = k, n = n, lambda = lambda_hat, theta = theta_hat, alpha = alpha_hat)
library(vegan)
specpool(fungi2016_mat)

rarefaction_pois(n, lambda = lambda_hat, theta = theta_hat, alpha = alpha_hat)
pois_fit <- rarefaction_pois(n, lambda = lambda_hat, theta = theta_hat, alpha = alpha_hat)
pois_fit

# IBP case
fit_IBP <- EB_EFPF_gamma_IBP(n, freq, pars_0 = c(0.1, 10, 200, 1))
alpha_hat <- fit_IBP[1]
theta_hat <- fit_IBP[2]
gamma_hat <- fit_IBP[3] / fit_IBP[4]
-alpha_hat / theta_hat
theta_hat
gamma_hat

IBP_fit <- rarefaction_IBP(n, gamma = gamma_hat, theta = theta_hat, alpha = alpha_hat)
IBP_fit


# ------------------------------

library(vegan)

# Accumulation curve using vegan
vegan_fit <- poolaccum(fungi2016_mat)

vegan_fit

# Plot of the raw date (random order + rarefaction)
plot_trajectory(fungi2016_list) + xlab("Plot") + ylab("Species") + ggtitle("") + 
  geom_line(data = data.frame(sites = c(0,vegan_fit$means[, 1]), richness = c(0,vegan_fit$means[, 2])), aes(x = sites, y = richness), linetype = "dotted") +
  geom_line(data = data.frame(sites = c(0,pois_fit$sites), richness = c(0,pois_fit$richness)), aes(x = sites, y = richness), linetype = "dashed", col = "darkblue")#+
#  geom_line(data = data.frame(sites = c(0, IBP_fit$sites), richness = c(0, IBP_fit$richness)), aes(x = sites, y = richness), linetype = "dashed", col = "darkorange")



