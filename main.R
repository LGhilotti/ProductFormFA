######################################################
#### Main script where functions are invoked ##########
#######################################################

source("R/routines_poiss_BB.R")
source("R/routines_negbin_BB.R")
source("R/routines_gamma_IBP.R")

## BB with Poisson(lambda)
set.seed(1234)
alpha <- -1
theta <- 15
lambda <- 100
lev <- 0.95
n <- 10
m <- 10

plot_Kmn_poiss_BB(alpha,theta,m,n,lambda, lev)

buff_poiss_bb <- buffet_poiss_BB(alpha,theta,n,lambda)

## BB with Negative-Binomial(n*,p) 
set.seed(1234)
alpha <- -1
theta <- 15
nstar <- 100
p <- 0.5
lev <- 0.95
n <- 10
Kn <- 5
m <- 10

plot_Kmn_negbin_BB(alpha,theta,m,n,Kn,nstar, p, lev)

buff_negbin_bb <- buffet_negbin_BB(alpha,theta,n,nstar,p)

## IBP with Gamma(a,b)
set.seed(1234)
alpha <- 0.2
theta <- 5
a <- 10
b <- 3
lev <- 0.95
n <- 10
Kn <- 3
m <- 10

plot_Kmn_gamma_IBP(alpha,theta,m,n,Kn,a, b, lev)

buff_gamma_ibp <- buffet_gamma_IBP(alpha,theta,n,a,b)
