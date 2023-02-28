######################################################
#### Main script where functions are invoked ##########
#######################################################

source("R/routines_poiss_BB.R")
source("R/routines_negbin_BB.R")

## BB with Poisson(lambda)
alpha <- -1
theta <- 15
lambda <- 100
lev <- 0.95
n <- 10
m <- 10

plot_Kmn_poiss_BB(alpha,theta,m,n,lambda, lev)

## BB with NegBin(n*,p) 
alpha <- -1
theta <- 15
nstar <- 10
p <- 0.5
lev <- 0.95
n <- 10
Kn <- 1
m <- 10

plot_Kmn_negbin_BB(alpha,theta,m,n,Kn,nstar, p, lev)


## IBP with Gamma(a,b)
