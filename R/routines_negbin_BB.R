#' Credible intervals of Kmn for the BB with Negative-Binomial mixture
#'
#' This function computes the means and the credible intervals of Kmn for the BB with Negative-Binomial mixture
#'
#' @param alpha [numeric] value of alpha in product-form feature allocation
#' @param theta [numeric] value of theta in product-form feature allocation
#' @param m [numeric] dimension of the new sample to be observed
#' @param n [numeric] dimension of the already observed sample
#' @param Kn [numeric] number of features in the already observed sample
#' @param nstar [numeric] Negative-Binomial hyperparameter (number of successes)
#' @param p [numeric] Negative-Binomial hyperparameter (success probability)
#' @param lev [numeric] level of the credible intervals
#'
#' @export
#'
CI_Kmn_negbin_BB <- function(alpha, theta, m, n, Kn, nstar, p, lev) {
  
  pbars <- p_kmn_all_negbin_BB(alpha, theta, m, n, p)

  means <- (nstar + Kn) * (1 - pbars) / pbars

  ub <- qnbinom(lev + (1 - lev) / 2, nstar + Kn, pbars)

  lb <- qnbinom((1 - lev) / 2, nstar + Kn, pbars)

  return(list("means" = means,"ubs" = ub,"lbs" = lb))
  
}


#########################################################################

#' EB based on EFPF-max (possibly some parameters are fixed) - 
#' - BB with Negative-Binomial mixture
#'
#' This function returns the value of the parameters maximizing the 
#' EFPF for the given sample (possibly some parameters are fixed) - 
#' - BB with Negative-Binomial mixture
#'
#' @param n [integer] dimension of the observed sample
#' @param counts [numeric] vector of cardinalities for the observed features
#' @param pars_0 [numeric] Initialization for (alpha, theta, n*, p) to 
#' be optimized
#' @param opt [logical] vector with 1 if value is optimized, 0 if not optimized
#'  
#' @export
EB_EFPF_fixed_negbin_BB <- function(n, counts, pars_0, opt = c(T,T,T,T)){
  
  lb <- c(-Inf,0.01, 0.01, 0.01)
  ub <- c(-0.01, Inf, Inf, 0.99)
  
  # if alpha or theta are fixed -> optimize the others: not need to reparametrize
  if (opt[1] ==F | opt[2] == F){
    if (opt[1] == F){
      lb[2] <- -pars_0[1]
    }
    else {
      lb[1] <- -pars_0[2]
    }
    res <- optim(par = pars_0[opt], fn = neg_log_EFPF_negbin_BB, n = n, 
                 counts = counts, nopt_par = pars_0[!opt], opt = opt,
                 method = "L-BFGS-B", lower = lb[opt], 
                 upper = ub[opt])
    sol <- pars_0
    sol[opt] <- res$par
  }
  else{
    pars_0[2] <- pars_0[2] + pars_0[1]
    # set constraints stricter so that function is limited
    if (sum(opt)==4){
      res <- optim(par = pars_0, fn = neg_log_EFPF_negbin_BB_rep_all, n = n, 
                   counts = counts,
                   method = "L-BFGS-B", lower = lb, 
                   upper = ub)
    }
    else{
      res <- optim(par = pars_0[opt], fn = neg_log_EFPF_negbin_BB_rep, n = n, 
                   counts = counts, nopt_par = pars_0[!opt], opt = opt,
                   method = "L-BFGS-B", lower = lb[opt], 
                   upper = ub[opt])
    }
    
    
    sol <- pars_0
    sol[opt] <- res$par
    sol[2] <- sol[2] - sol[1]
    
  }
    
  
  return (sol)
  
}


#########################################################
###### mean-variance NegBin parametrization ############
########################################################
#########################################################################

neg_log_EFPF_negbin_mv_BB_rep_R <- function(pars, n, counts, nopt_par, opt){
  
  i_T=1; i_F=1;
  # Assign alpha
  if (opt[1]==TRUE){
    alpha = pars[1];
    i_T = i_T +1;
    # Assign s
    if (opt[2]==TRUE){
      theta = pars[i_T] - alpha;
      i_T = i_T +1;
    }
    else{
      theta = nopt_par[i_F];
      i_F = i_F +1;
    }
  }
  else {
    alpha = nopt_par[1]; 
    i_F = i_F +1;
    if (opt[2]==TRUE){
      theta = pars[i_T];
      i_T = i_T +1;
    }
    else{
      theta = nopt_par[i_F];
      i_F = i_F +1;
    }
  }
  
  # Assign mu
  if (opt[3]==TRUE){
    mu = pars[i_T];
    i_T = i_T +1;
    # Assign t
    if (opt[4]==TRUE){
      sigma2 = pars[i_T] + mu;
    }
    else{
      sigma2 = nopt_par[i_F];
    }
  }
  else {
    mu = nopt_par[i_F]; 
    i_F = i_F +1;
    # Assign t
    if (opt[4]==TRUE){
      sigma2 = pars[i_T];
    }
    else{
      sigma2 = nopt_par[i_F];
    }
  }
  
  s <- alpha + theta
  t <- sigma2 - mu
  
  return (neg_log_EFPF_negbin_mv_BB_rep_all(pars = c(alpha,s,mu,t), n, counts))
}

#' EB based on EFPF-max (possibly some parameters are fixed) - 
#' - BB with Negative-Binomial mixture (mean-variance parametrization)
#'
#' This function returns the value of the parameters maximizing the 
#' EFPF for the given sample (possibly some parameters are fixed) - 
#' - BB with Negative-Binomial mixture (mean-variance parametrization)
#'
#' @param n [integer] dimension of the observed sample
#' @param counts [numeric] vector of cardinalities for the observed features
#' @param pars_0 [numeric] Initialization for (alpha, theta, mu, sigma2) to 
#' be optimized
#' @param opt [logical] vector with 1 if value is optimized, 0 if not optimized
#'  
#' @export
EB_EFPF_fixed_negbin_mv_BB <- function(n, counts, pars_0, opt = c(T,T,T,T)){
  
  lb <- c(-Inf,0.1, 0.1, 0.1)
  ub <- c(-0.1, Inf, Inf, Inf)
  
  # if alpha or theta are fixed AND mu or sigma2 are fixed -> 
  # -> optimize the others: not need to reparametrize
  if ((opt[1] ==F | opt[2] == F) & (opt[3] ==F | opt[4] == F) ){
    if (opt[1] == F){
      lb[2] <- -pars_0[1]+0.01
    }
    else {
      lb[1] <- -pars_0[2]+0.01
    }
    
    if (opt[3] == F){
      lb[4] <- pars_0[3]+0.01
    }
    else{
      ub[3] <- pars_0[4]-0.01
    }
    res <- optim(par = pars_0[opt], fn = neg_log_EFPF_negbin_mv_BB, n = n, 
                 counts = counts, nopt_par = pars_0[!opt], opt = opt,
                 method = "L-BFGS-B", lower = lb[opt], 
                 upper = ub[opt])
    sol <- pars_0
    sol[opt] <- res$par
  }
  else if (sum(opt)==4){
    pars_0[2] <- pars_0[1] + pars_0[2]
    pars_0[4] <- pars_0[4] - pars_0[3]
    # set constraints stricter so that function is limited
    res <- optim(par = pars_0, fn = neg_log_EFPF_negbin_mv_BB_rep_all, n = n, 
                 counts = counts,
                 method = "L-BFGS-B", lower = lb, upper = ub)
    
    sol <- pars_0
    sol[opt] <- res$par
    sol[2] <- sol[2] - sol[1]
    sol[4] <- sol[4] + sol[3]
  }
  else {
    if (opt[1] == T & opt[2] == T){ #alpha or theta are optimized -> use s
      pars_0[2] <- pars_0[1] + pars_0[2]
    }
    if (opt[3] == T & opt[4] == T){ #mu or sigma are optimized -> use t
      pars_0[4] <- pars_0[4] - pars_0[3]
    }  
    if (opt[3] == F){
      lb[4] <- pars_0[3]+0.01
    }
    if (opt[4] == F){
      ub[3] <- pars_0[4]-0.01
    }
      
    if (opt[1] == F){
      lb[2] <- -pars_0[1]+0.01
    }
    if (opt[2] == F) {
      lb[1] <- -pars_0[2]+0.01
    }
    res <- optim(par = pars_0[opt], fn = neg_log_EFPF_negbin_mv_BB_rep_R, n = n, 
                 counts = counts, nopt_par = pars_0[!opt], opt = opt,
                 method = "L-BFGS-B", lower = lb[opt], 
                 upper = ub[opt])
    
    sol <- pars_0
    sol[opt] <- res$par
    if (opt[1] == T & opt[2] == T){ #alpha or theta are optimized -> use s
      sol[2] <- sol[2] - sol[1]
    }
    if (opt[3] == T & opt[4] == T){ #mu or sigma are optimized -> use t
      sol[4] <- sol[4] + sol[3]
    }
  }
  
  
  return (sol)
  
}
