
########### BB with Poisson(lambda) mixture ###############

neg_log_posterior_poiss <- function(pars,
                                    n, K, counts,
                                    lambda,
                                    a_s, b_s, a_alpha, b_alpha){
  
  s_hat <- pars[1]
  alpha_hat <- pars[2]
  
  es <- exp(s_hat)
  ealpha <- exp(alpha_hat)
  
  res <- K*(lgamma(es +ealpha) -lgamma(es + ealpha +n) -
              lgamma(1 + ealpha) -lgamma(es) ) +
    lambda*exp( lgamma(es+n) +lgamma(es+ealpha) - lgamma(es) - lgamma(es+ealpha+n))  +
    sum(lgamma(ealpha + counts) +lgamma(es+n-counts)) +
    alpha_hat*(K+a_alpha) - b_alpha*ealpha +
    s_hat*a_s - b_s*es
  
  return (-res)
}


compute_grad_log_full_poiss <- function(s_hat, alpha_bar_hat, lambda,
                                        n, K, counts, a_s, b_s, a_alpha, b_alpha){
  
  es <- exp(s_hat)
  ea_bar <- exp(alpha_bar_hat)
  
  p1 <- exp(lgamma(n+es) + lgamma(es + ea_bar) -
              lgamma(es) - lgamma(es + ea_bar +n) )
  
  ds_hat <- a_s + es * (lambda*p1*(digamma(n+es) - digamma(es) - digamma(n + es + ea_bar) + 
                                     digamma(es + ea_bar)) - 
                          K*(digamma(n+es+ea_bar) - digamma(es+ea_bar)) +
                          sum(digamma(n - counts + es)) -
                          K*digamma(es) - b_s)
  
  dalpha_bar_hat <- K + a_alpha + ea_bar*( (digamma(n+es+ea_bar) - digamma(es+ea_bar))*
                                             (- lambda* p1 - K) +
                                             sum(digamma(counts + ea_bar)) -
                                             K*digamma(1+ea_bar) - b_alpha)
  
  return (c(ds_hat, dalpha_bar_hat))
}


compute_log_ratio_q_poiss_precond <- function(params_prop, params_curr,
                                              tau, cov_post,
                                              grad_log_full_curr, grad_log_full_prop){
  
  norm_prop <- t(params_curr-params_prop-tau*cov_post%*%grad_log_full_prop)%*%
    inv(cov_post)%*%(params_curr-params_prop-tau*cov_post%*%grad_log_full_prop)
  
  norm_curr <- t(params_prop-params_curr-tau*cov_post%*%grad_log_full_curr)%*%
    inv(cov_post)%*%(params_prop-params_curr-tau*cov_post%*%grad_log_full_curr)
  
  res <- 1/(4*tau)*(norm_curr - norm_prop)
  
  return (res)
}


#' Metropolis-within-Gibbs sampler for BB with Poisson(lambda) mixture
#'
#' @param Z [integer] binary matrix of presence/absence (n x K - dimensional)
#' @param alpha_bar_0 [numeric] initial value of alpha_bar 
#' @param s_0 [numeric] initial value of s
#' @param a_alpha [numeric]
#' @param b_alpha [numeric]
#' @param a_s [numeric]
#' @param b_s [numeric]
#' @param lambda [numeric] 
#' @param tau [numeric] MALA step-size
#' @param S [integer] number of iterations for the MCMC algorithm
#' @param n_burnin [integer] number of iterations for the burn-in
#' @param thin [integer] thinning
#' @param seed [integer] seed
#'
#' @return
#'
#' @import numDeriv
#' @import stats
#' @import MASS
#' @import matlib
#' 
sampler_PoissonBB <- function(Z,
                              alpha_bar_0, s_0,
                              a_alpha, b_alpha, a_s, b_s, lambda,
                              tau, S, n_burnin, thin, seed){
  
  set.seed(seed)
  
  # Compute total number of sites
  n <- nrow(Z)
  
  # Delete NA
  Z <- Z[, colSums(is.na(Z))==0]
  
  # Delete zero-columns
  Z <- Z[, colSums(Z)!=0]
  
  # Set K to be the observed number of features
  K <- ncol(Z)
  
  # Compute vector of counts
  counts <- colSums(Z)
  
  ############## Gibbs-sampler ##########################
  
  # Define structure to store parameters along the iterations
  number_saved_iterations <- (S - n_burnin)/thin + 1
  alpha_bar_vec <- vector(length = number_saved_iterations)
  s_vec <- vector(length = number_saved_iterations)
  
  # Set initial values
  alpha_bar <- alpha_bar_0
  s <- s_0
  
  # Compute mode of the log-posterior density (order: s_hat, alpha_hat)
  mode_post <- optim(par = c(0,0), fn = neg_log_posterior_poiss, n = n, K = K, counts = counts,
                     lambda = lambda, a_s = a_s, b_s = b_s, a_alpha = a_alpha, b_alpha = b_alpha,
                     method = "L-BFGS-B")$par
  
  # Compute the hessian of the neg log-posterior in the mode
  hess_neg_log_mode <- hessian(func = neg_log_posterior_poiss, x=mode_post, 
                               n = n, K = K, counts = counts,
                               lambda = lambda, a_s = a_s, b_s = b_s, a_alpha = a_alpha, b_alpha = b_alpha)
  
  # Covariance matrix of posterior density
  cov_post <- inv(hess_neg_log_mode)
  
  # index saved iterations (after burn-in and thinning satisfied)
  l <- 1
  
  for (q in 1:S){
    
    
    ################################################################
    ############# Draw (s, alpha_bar) | Z ##################
    ###############################################################
    
    # In order to update s, alpha_bar, we update s_hat, alpha_bar_hat, defined
    # as the logarithm of s, alpha_bar
    
    ### Current values for s_hat, alpha_bar_hat
    s_hat_curr <- log(s)
    alpha_bar_hat_curr <- log(alpha_bar)
    
    ### Propose values for s_hat, alpha_bar_hat
    # Compute the gradient of log-full conditional for MALA
    grad_log_full_curr <- compute_grad_log_full_poiss(s_hat_curr, alpha_bar_hat_curr, lambda,
                                                      n, K, counts, a_s, b_s, a_alpha, b_alpha)
    
    # s_hat_prop <- s_hat_curr + tau*grad_log_full_curr[1] + sqrt(2*tau)*rnorm(1)
    # alpha_bar_hat_prop <- alpha_bar_hat_curr + tau*grad_log_full_curr[2] + sqrt(2*tau)*rnorm(1)
    
    # Propose from the bivariate normal 
    params_curr <- c(s_hat_curr, alpha_bar_hat_curr)
    params_prop <- mvrnorm(mu = params_curr + tau*cov_post%*%grad_log_full_curr,
                           Sigma = 2*tau*cov_post)
    
    s_hat_prop <- params_prop[1]
    alpha_bar_hat_prop <- params_prop[2]
    
    ### Acceptance probability 
    # Compute the log ratio of the full-cond in prop point and curr point
    # log_ratio_full <- compute_log_ratio_full_poiss(s_hat_prop, alpha_bar_hat_prop,
    #                                                s_hat_curr, alpha_bar_hat_curr, lambda,
    #                                                n, K, counts, a_s, b_s, a_alpha, b_alpha)
    log_ratio_full <- - neg_log_posterior_poiss(params_prop, n, K, counts,
                                                lambda, a_s, b_s, a_alpha, b_alpha) +
      neg_log_posterior_poiss(params_curr, n, K, counts,
                              lambda, a_s, b_s, a_alpha, b_alpha )
    
    # Compute the log ratio of the terms related to the proposal q
    grad_log_full_prop <- compute_grad_log_full_poiss(s_hat_prop, alpha_bar_hat_prop, lambda,
                                                      n, K, counts, a_s, b_s, a_alpha, b_alpha)
    
    # log_ratio_q <- compute_log_ratio_q_poiss(s_hat_prop, alpha_bar_hat_prop,
    #                                          s_hat_curr, alpha_bar_hat_curr, tau,
    #                                          grad_log_full_curr, grad_log_full_prop)
    log_ratio_q <- compute_log_ratio_q_poiss_precond(params_prop, params_curr,
                                                     tau, cov_post,
                                                     grad_log_full_curr, grad_log_full_prop)
    
    # Compute acceptance probability
    log_acc_prob <- log_ratio_full + log_ratio_q
    acc_prob <- min(1, exp(log_acc_prob))
    
    # Decide if accept or not the new parameter vector
    if (runif(1) < acc_prob){ # accept
      s <- exp(s_hat_prop)
      alpha_bar <- exp(alpha_bar_hat_prop)
    }
    
    
    #######################################################################
    
    # Store parameters if burn-in is over and once every "thin" iteration
    if ((q > n_burnin) & (q %% thin == 0) ){
      print(paste0("iteration: ", q))
      
      alpha_bar_vec[l] <- alpha_bar
      s_vec[l] <- s
      
      l <- l+1
    }
    
  }
  
  alpha_bar_vec <- alpha_bar_vec[1:(l-1)]
  s_vec <- s_vec[1:(l-1)]
  
  return (list("alpha_bar_chain" = alpha_bar_vec, 
               "s_chain" = s_vec))
  
}









##################### BB with NB(n0,mu0) mixture #######################

neg_log_posterior_negbin <- function(pars,
                                     n, K, counts,
                                     nstar, p,
                                     a_s, b_s, a_alpha, b_alpha){
  
  s_hat <- pars[1]
  alpha_hat <- pars[2]
  
  es <- exp(s_hat)
  ealpha <- exp(alpha_hat)
  
  res <- K*(lgamma(es +ealpha) -lgamma(es + ealpha +n) -
              lgamma(1 + ealpha) -lgamma(es) ) -
    (K+nstar)*log(1-(1-p)*exp( lgamma(es+n) +lgamma(es+ealpha) - lgamma(es) - lgamma(es+ealpha+n)) ) +
    sum(lgamma(ealpha + counts) +lgamma(es+n-counts)) +
    alpha_hat*(K+a_alpha) - b_alpha*ealpha +
    s_hat*a_s - b_s*es
  
  return (-res)
}


compute_grad_log_full_negbin <- function(s_hat, alpha_bar_hat, nstar, p,
                                              n, K, counts, 
                                              a_s, b_s, a_alpha, b_alpha){
  
  es <- exp(s_hat)
  ea_bar <- exp(alpha_bar_hat)
  
  p1 <- exp(lgamma(n+es) + lgamma(es + ea_bar) -
              lgamma(es) - lgamma(es + ea_bar +n) )
  
  # derivative wrt s_hat
  
  p2 <- exp( log(nstar + K) + log(1-p) - log(1/p1 -1 + p) )
  
  ds_hat <- a_s + es * (p2*(digamma(n+es) - digamma(es) - digamma(n + es + ea_bar) + 
                              digamma(es + ea_bar)) - 
                          K*(digamma(n+es+ea_bar) - digamma(es+ea_bar)) +
                          sum(digamma(n - counts + es)) -
                          K*digamma(es) - b_s)
  
  # derivative wrt alpha_bar_hat
  dalpha_bar_hat <- K + a_alpha + ea_bar*( (digamma(n+es+ea_bar) - digamma(es+ea_bar))*
                                             (- p2 - K) +
                                             sum(digamma(counts + ea_bar)) -
                                             K*digamma(1+ea_bar) - b_alpha)
  
  return (c(ds_hat, dalpha_bar_hat))
}

compute_log_ratio_q_negbin_precond <- function(params_prop, params_curr,
                                               tau, cov_post,
                                               grad_log_full_curr, grad_log_full_prop){
  
  norm_prop <- t(params_curr-params_prop-tau*cov_post%*%grad_log_full_prop)%*%
    inv(cov_post)%*%(params_curr-params_prop-tau*cov_post%*%grad_log_full_prop)
  
  norm_curr <- t(params_prop-params_curr-tau*cov_post%*%grad_log_full_curr)%*%
    inv(cov_post)%*%(params_prop-params_curr-tau*cov_post%*%grad_log_full_curr)
  
  res <- 1/(4*tau)*(norm_curr - norm_prop)
  
  return (res)
}

#' Metropolis-within-Gibbs sampler for BB with NB(n0,mu0) mixture
#'
#' @param Z [integer] binary matrix of presence/absence (n x K- dimensional)
#' @param alpha_bar_0 [numeric] initial value of alpha_bar 
#' @param s_0 [numeric] initial value of s
#' @param a_alpha [numeric]
#' @param b_alpha [numeric]
#' @param a_s [numeric]
#' @param b_s [numeric]
#' @param nstar [numeric] 
#' @param p [numeric] 
#' @param tau [numeric] MALA step-size
#' @param S [integer] number of iterations for the MCMC algorithm
#' @param n_burnin [integer] number of iterations for the burn-in
#' @param thin [integer] thinning
#' @param seed [integer] seed
#'
#' @return
#' 
#' @import numDeriv
#' @import stats
#' @import MASS
#' @import matlib
#'
sampler_NegBinBB <- function(Z, 
                             alpha_bar_0, s_0,
                             a_alpha, b_alpha, a_s, b_s, nstar, p,
                             tau, S, n_burnin, thin, seed){
  
  set.seed(seed)
  
  # Compute total number of sites
  n <- nrow(Z)
  
  # Delete NA
  Z <- Z[, colSums(is.na(Z))==0]
  
  # Delete zero-columns
  Z <- Z[, colSums(Z)!=0]
  
  # Set K to be the observed number of features
  K <- ncol(Z)
  
  # Compute vector of counts
  counts <- colSums(Z)
  
  ############## Gibbs-sampler ##########################
  
  # Define structure to store parameters along the iterations
  number_saved_iterations <- (S - n_burnin)/thin + 1
  alpha_bar_vec <- vector(length = number_saved_iterations)
  s_vec <- vector(length = number_saved_iterations)
  
  # Set initial values
  alpha_bar <- alpha_bar_0
  s <- s_0
  
  # Compute mode of the log-posterior density (order: s_hat, alpha_hat)
  mode_post <- optim(par = c(0,0), fn = neg_log_posterior_negbin, n = n, K = K, counts = counts,
                     nstar = nstar, p=p,
                     a_s = a_s, b_s = b_s, a_alpha = a_alpha, b_alpha = b_alpha,
                     method = "L-BFGS-B")$par
  
  # Compute the hessian of the neg log-posterior in the mode
  hess_neg_log_mode <- hessian(func = neg_log_posterior_negbin, x=mode_post, 
                               n = n, K = K, counts = counts,
                               nstar = nstar, p=p, a_s = a_s, b_s = b_s, a_alpha = a_alpha, b_alpha = b_alpha)
  
  # Covariance matrix of posterior density
  cov_post <- inv(hess_neg_log_mode)
  
  # index saved iterations (after burn-in and thinning satisfied)
  l <- 1
  
  for (q in 1:S){
    
    t <- lgamma(s + n) + lgamma(s + alpha_bar) - lgamma(s) - lgamma(s+alpha_bar+n)
    
    ################################################################
    ############# Draw ( s, alpha_bar) | Z  ##################
    ###############################################################
    
    # In order to update s, alpha_bar, we update s_hat, alpha_bar_hat, 
    # defined as the logarithm of s, alpha_bar
    
    ### Current values for s_hat, alpha_bar_hat
    s_hat_curr <- log(s)
    alpha_bar_hat_curr <- log(alpha_bar)
    
    ### Propose values for s_hat, alpha_bar_hat
    # Compute the gradient of log-full conditional for MALA 
    # (the order of returned variables is: s_hat, alpha_bar_hat)
    grad_log_full_curr <- compute_grad_log_full_negbin(s_hat_curr, alpha_bar_hat_curr, 
                                                            nstar, p, n, K, counts, 
                                                            a_s, b_s, a_alpha, b_alpha)
    
    # s_hat_prop <- s_hat_curr + tau*grad_log_full_curr[1] + sqrt(2*tau)*rnorm(1)
    # 
    # alpha_bar_hat_prop <- alpha_bar_hat_curr + tau*grad_log_full_curr[2] + sqrt(2*tau)*rnorm(1)
    
    # Propose from the bivariate normal 
    params_curr <- c(s_hat_curr, alpha_bar_hat_curr)
    params_prop <- mvrnorm(mu = params_curr + tau*cov_post%*%grad_log_full_curr,
                           Sigma = 2*tau*cov_post)
    
    s_hat_prop <- params_prop[1]
    alpha_bar_hat_prop <- params_prop[2]
    
    ### Acceptance probability 
    # Compute the log ratio of the full-cond in prop point and curr point
    # log_ratio_full <- compute_log_ratio_full_negbin_geom(s_hat_prop, alpha_bar_hat_prop,
    #                                                      s_hat_curr, alpha_bar_hat_curr, 
    #                                                      nstar, p, n, K, counts,   
    #                                                      a_s, b_s, a_alpha, b_alpha)
    log_ratio_full <- - neg_log_posterior_negbin(params_prop, n, K, counts,
                                                 nstar,p, a_s, b_s, a_alpha, b_alpha) +
      neg_log_posterior_negbin(params_curr, n, K, counts,
                               nstar, p, a_s, b_s, a_alpha, b_alpha )
    
    # Compute the log ratio of the terms related to the proposal q
    grad_log_full_prop <- compute_grad_log_full_negbin(s_hat_prop, alpha_bar_hat_prop, 
                                                            nstar, p, n, K, counts, 
                                                            a_s, b_s, a_alpha, b_alpha)
    
    # log_ratio_q <- compute_log_ratio_q_negbin_geom(s_hat_prop, alpha_bar_hat_prop,
    #                                                s_hat_curr, alpha_bar_hat_curr, tau,
    #                                                grad_log_full_curr, grad_log_full_prop)
    log_ratio_q <- compute_log_ratio_q_negbin_precond(params_prop, params_curr,
                                                      tau, cov_post,
                                                      grad_log_full_curr, grad_log_full_prop)
    
    # Compute acceptance probability
    log_acc_prob <- log_ratio_full + log_ratio_q
    acc_prob <- min(1, exp(log_acc_prob))
    
    # Decide if accept or not the new parameter vector
    if (runif(1) < acc_prob){ # accept
      s <- exp(s_hat_prop)
      alpha_bar <- exp(alpha_bar_hat_prop)
    }
    
    
    #######################################################################
    
    # Store parameters if burn-in is over and once every "thin" iteration
    if ((q > n_burnin) & (q %% thin == 0) ){
      print(paste0("iteration: ", q))
      
      alpha_bar_vec[l] <- alpha_bar
      s_vec[l] <- s
      
      l <- l+1
    }
    
  }
  
  alpha_bar_vec <- alpha_bar_vec[1:(l-1)]
  s_vec <- s_vec[1:(l-1)]
  
  return (list("alpha_bar_chain" = alpha_bar_vec, "s_chain" = s_vec))
  
}




################ IBP with Gamma(a,b) mixture ####################

neg_log_posterior_gamma_ibp <- function(pars, 
                                        n, K, counts,
                                        a,b,
                                        a_s, b_s, a_alpha, b_alpha){
  
  s_hat <- pars[1]
  alpha_hat <- pars[2]
  
  es <- exp(s_hat)
  ealpha <- exp(alpha_hat)
  
  res <- K*(lgamma(es - ealpha/(1+ealpha) +1) -lgamma(es - ealpha/(1+ealpha) +n) -
              lgamma(1 - ealpha/(1+ealpha)) -lgamma(es) ) -
    (K+a)*log(b + exp(lgamma(es - ealpha/(1+ealpha) +1) -lgamma(es))*
                sum(exp(lgamma(es + 1:n -1) -lgamma(es - ealpha/(1+ealpha) +1:n))) ) +
    sum(lgamma(-ealpha/(1+ealpha) + counts) +lgamma(es+n-counts)) +
    alpha_hat*a_alpha - (a_alpha + b_alpha)*log(1+ealpha) +
    s_hat*a_s - b_s*es
  
  return (-res)
}

#' Metropolis-within-Gibbs sampler for IBP with Gamma(a,b) mixture,
#' with prior also on a and b
#'
#' @param Z [integer] binary matrix of presence/absence (n x K- dimensional)
#' @param alpha_0 [numeric] initial value of alpha 
#' @param s_0 [numeric] initial value of s
#' @param a_0 [numeric] initial value of a
#' @param b_0 [numeric] initial value of b
#' @param a_alpha [numeric]
#' @param b_alpha [numeric]
#' @param a_s [numeric]
#' @param b_s [numeric]
#' @param q [numeric] hyperparameter of Geometric prior on a
#' @param r [numeric] 
#' @param t [numeric]
#' @param tau [numeric] MALA step-size
#' @param fixed [logical] vector indicating if the variable is fixed or has prior:
#' if fixed, it stays equal to initial value - order: a, b, alpha, s
#' @param S [integer] number of iterations for the MCMC algorithm
#' @param n_burnin [integer] number of iterations for the burn-in
#' @param thin [integer] thinning
#' @param seed [integer] seed
#'
#' @return
#' 
#' @import numDeriv
#' @import stats
#' @import MASS
#' @import matlib
#'
sampler_GammaIBP <- function(Z, 
                             alpha_0, s_0, a_0, b_0,
                             a_alpha, b_alpha, a_s, b_s, q, r, t,
                             sigq_alpha, sigq_s, S, n_burnin, thin, seed){
  
  set.seed(seed)
  
  # Compute total number of sites
  n <- nrow(Z)
  
  # Delete NA
  Z <- Z[, colSums(is.na(Z))==0]
  
  # Delete zero-columns
  Z <- Z[, colSums(Z)!=0]
  
  # Set K to be the observed number of features
  K <- ncol(Z)
  
  # Compute vector of counts
  counts <- colSums(Z)
  
  ############## Gibbs-sampler ##########################
  
  # Define structure to store parameters along the iterations
  number_saved_iterations <- (S - n_burnin)/thin + 1
  a_vec <- vector(length = number_saved_iterations)
  b_vec <- vector(length = number_saved_iterations)
  alpha_vec <- vector(length = number_saved_iterations)
  s_vec <- vector(length = number_saved_iterations)
  
  # Set initial values
  a <- a_0
  b <- b_0
  alpha <- alpha_0
  s <- s_0
  
  # index saved iterations (after burn-in and thinning satisfied)
  l <- 1
  
  for (q in 1:S){
    
    ################################################################
    ############# Draw a, b, s, alpha | Z   ###################
    ###############################################################
    
    # Update Gam | Z, alpha, s, a, b
    
    gamma_a_s_n <- sum(exp(lgamma(s + (1:n) - 1) - lgamma(s) - 
                             lgamma(s - alpha + (1:n)) + lgamma(s - alpha +1) ) )
    
    Gam <- rgamma(1, shape = K + a, rate = b + gamma_a_s_n)
    
   
    # Update a | Gam, Z, alpha, s, b
    
    a <- 1 + rpois(1, b*Gam*(1-q) )
    

    # Update b | Gam, Z, alpha, s, a
    
    b <- rgamma(1, shape = a + r, rate = Gam + t )
    
    
    ################################################################
    ############# Draw ( s, alpha) | Z, a, b ##################
    ###############################################################
    
    # In order to update s, alpha, we update s_hat, alpha_hat, 
    # defined as the logarithm of s, alpha
    
    ### Current values for s_hat, alpha_hat
    s_hat_curr <- log(s)
    alpha_hat_curr <- log(alpha/(1-alpha))
    params_curr <- c(s_hat_curr, alpha_hat_curr)
    
    ### Update s_hat first
    # Propose values for s_hat 
    s_hat_prop <- rnorm(n=1, mean = s_hat_curr, sd = sqrt(sigq_s) ) 
    
    params_prop <- c(s_hat_prop, alpha_hat_curr)
    
    # Compute acceptance probability: log ratio of the full-cond in prop point and curr point
    log_acc_prob <- - neg_log_posterior_gamma_ibp(params_prop, n, K, counts,
                                                  a,b, a_s, b_s, a_alpha, b_alpha) +
      neg_log_posterior_gamma_ibp(params_curr, n, K, counts,
                                  a,b, a_s, b_s, a_alpha, b_alpha )
    
    acc_prob <- min(1, exp(log_acc_prob))
    
    # Decide if accept or not the new parameter 
    if (runif(1) < acc_prob){ # accept
      s_hat_curr <- s_hat_prop
      params_curr <- params_prop
    }
    
    ### Update alpha_hat then
    alpha_hat_prop <- rnorm(n=1, mean = alpha_hat_curr, sd = sqrt(sigq_alpha) )
   
    params_prop <- c(s_hat_curr, alpha_hat_prop)
    
    # Compute acceptance probability: log ratio of the full-cond in prop point and curr point
    log_acc_prob <- - neg_log_posterior_gamma_ibp(params_prop, n, K, counts,
                                                  a,b, a_s, b_s, a_alpha, b_alpha) +
      neg_log_posterior_gamma_ibp(params_curr, n, K, counts,
                                  a,b, a_s, b_s, a_alpha, b_alpha )
    
    acc_prob <- min(1, exp(log_acc_prob))
    
    # Decide if accept or not the new parameter 
    if (runif(1) < acc_prob){ # accept
      alpha_hat_curr <- alpha_hat_prop
    }
    
    
    ### Set the vector (s, alpha) with the updated values
    s <- exp(s_hat_curr)
    alpha <- exp(alpha_hat_curr)/(1 + exp(alpha_hat_curr))
    
    
    #######################################################################
    
    # Store parameters if burn-in is over and once every "thin" iteration
    if ((q > n_burnin) & (q %% thin == 0) ){
      print(paste0("iteration: ", q))
      
      a_vec[l] <- a
      b_vec[l] <- b
      alpha_vec[l] <- alpha
      s_vec[l] <- s
      
      l <- l+1
    }
    
  }
  
  a_vec <- a_vec[1:(l-1)]
  b_vec <- b_vec[1:(l-1)]
  alpha_vec <- alpha_vec[1:(l-1)]
  s_vec <- s_vec[1:(l-1)]
  
  return (list("alpha_chain" = alpha_vec, "s_chain" = s_vec,
               "a_chain" = a_vec, "b_chain" = b_vec ))
  
}
