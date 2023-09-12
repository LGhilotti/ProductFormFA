rm(list=ls())
library(Rcpp)
library(matlib)
library(MASS)
library(numDeriv)

################ Utils function ##########################
create_features_list <- function(mat){
  
  feat_list <- vector("list", nrow(mat))
  
  for (i in 1:nrow(mat)){
    feat_list[[i]] <- which(mat[i,]==1, arr.ind = TRUE)
  }
  
  return (feat_list)
}

beta_binomial_estimator <- function(data_mat){
  
  # Compute total number of sites
  n <- nrow(data_mat)
  
  # Delete NA
  data_mat <- data_mat[, colSums(is.na(data_mat))==0]
  
  # Delete zero-columns
  data_mat <- data_mat[, colSums(data_mat)!=0]
  
  # Set K to be the observed number of features
  K <- ncol(data_mat)
  
  # Compute number of species contained in exactly k individuals
  counts <- colSums(data_mat)
  
  Q_1 <- sum(counts == 1)
  Q_2 <- sum(counts == 2)
  Q_3 <- sum(counts == 3)
  
  # Compute Q_hat_0 
  if (Q_2 == 0){
    Q_hat_0 <- (n-1)/n * (Q_1*(Q_1 -1))/2
  } else {
    Q_hat_0 <- (n-1)/n * (Q_1^2)/(2*Q_2)
  }
  
  # Compute last term
  if (2*Q_2^2 / (3*Q_1*Q_3) <= 1){
    last <- 2 - max(0.5, 2*Q_2^2 / (3*Q_1*Q_3))
  } else {
    last <- 1
  }
  
  # Compute the statistic
  res <- K + Q_hat_0 * last
  
  return (res)
  
}

#################### Functions for Poisson model #############################

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



compute_log_ratio_full_poiss <- function(s_hat_prop, alpha_bar_hat_prop,
                                         s_hat_curr, alpha_bar_hat_curr, lambda,
                                         n, K, counts, a_s, b_s, a_alpha, b_alpha){
  
  es_prop <- exp(s_hat_prop)
  ea_bar_prop <- exp(alpha_bar_hat_prop)
  es_curr <- exp(s_hat_curr)
  ea_bar_curr <- exp(alpha_bar_hat_curr)
  
  p1_prop <- exp(lgamma(n+es_prop) + lgamma(es_prop + ea_bar_prop) -
                   lgamma(es_prop) - lgamma(es_prop + ea_bar_prop +n) )
  p1_curr <- exp(lgamma(n+es_curr) + lgamma(es_curr + ea_bar_curr) -
                   lgamma(es_curr) - lgamma(es_curr + ea_bar_curr +n) )
  
  res <- lambda*(p1_prop - p1_curr) + (K+a_alpha)*(alpha_bar_hat_prop - alpha_bar_hat_curr) +
    a_s*(s_hat_prop - s_hat_curr) - b_alpha*(ea_bar_prop - ea_bar_curr) -
    b_s*(es_prop - es_curr) 
  
  # need to still add two terms
  res <- res + K*(lgamma(n+es_curr+ea_bar_curr) + lgamma(es_prop + ea_bar_prop) -
                    lgamma(es_curr + ea_bar_curr) - lgamma(es_prop+ea_bar_prop+n))
  
  v <- lgamma(ea_bar_prop + counts) - lgamma(1 + ea_bar_prop) + lgamma(es_prop+n - counts) -
    lgamma(es_prop) - lgamma(ea_bar_curr + counts) + lgamma(1 + ea_bar_curr) - 
    lgamma(es_curr+n - counts) + lgamma(es_curr)
  res <- res + sum(v)
  
  return(res)
  
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


compute_log_ratio_q_poiss <- function(s_hat_prop, alpha_bar_hat_prop,
                                      s_hat_curr, alpha_bar_hat_curr, tau,
                                      grad_log_full_curr, grad_log_full_prop){
  
  norm_prop <- norm(c(s_hat_curr, alpha_bar_hat_curr) - c(s_hat_prop, alpha_bar_hat_prop) -
                      tau*grad_log_full_prop, type="2")^2
  
  norm_curr <- norm( c(s_hat_prop, alpha_bar_hat_prop) - c(s_hat_curr, alpha_bar_hat_curr) -
                       tau*grad_log_full_curr, type="2")^2
  
  res <- 1/(4*tau)*(norm_curr - norm_prop)
  
  return (res)
}


gibbs_sampler_poiss_fixed_lambda <- function(Z, 
                                             alpha_bar_0, s_0,
                                             lambda, a_alpha, b_alpha, a_s, b_s,
                                             tau,
                                             S, n_burnin, thin, seed){
  
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
  
  return (list("alpha_bar_vec" = alpha_bar_vec, 
               "s_vec" = s_vec))
  
}




generate_Kmn_chain_poiss <- function(lambda, alpha_chain, theta_chain, M, n){
  
  S <- length(alpha_chain)
  
  if (length(lambda) == 1){
    lambda_chain <- rep(lambda, S)
  } else {
    lambda_chain <- lambda
  }
  
  M_vec <- 1:M
  kmn_chain <- matrix(NA, nrow = M, ncol = S )
  for (q in 1:S){
    lambda <- lambda_chain[q]
    alpha <- alpha_chain[q]
    theta <- theta_chain[q]
    
    par_1 <- lambda*exp(lgamma(theta+alpha+n) - lgamma(theta+alpha) - 
                          lgamma(theta+n) + lgamma(theta))
    par_2 <- lambda*exp(lgamma(theta+alpha+n+M_vec) - lgamma(theta+alpha+n) - 
                          lgamma(theta+n+M_vec) + lgamma(theta+n) +
                          lgamma(theta+alpha+n) - lgamma(theta+alpha) - 
                          lgamma(theta+n) + lgamma(theta) )
    poiss_par <- par_1 - par_2
    kmn_chain[,q] <- rpois(M, poiss_par)
  }
  
  return (kmn_chain)
}


generate_Ntilde_chain_poiss <- function(lambda, alpha_chain_poiss,
                                        theta_chain_poiss, n, Kn){
  S <- length(alpha_chain_poiss)
  
  if (length(lambda) == 1){
    lambda_chain_poiss <- rep(lambda, S)
  } else {
    lambda_chain_poiss <- lambda
  }  
  
  Ntilde_chain <- vector(length = S )
  for (q in 1:S){
    lambda <- lambda_chain_poiss[q]
    alpha <- alpha_chain_poiss[q]
    theta <- theta_chain_poiss[q]
    
    poiss_par <- lambda*exp(lgamma(theta+alpha+n) - lgamma(theta+alpha) - 
                              lgamma(theta+n) + lgamma(theta))
    
    Ntilde_chain[q] <- Kn + rpois(1, poiss_par)
  }
  
  return (Ntilde_chain)
}


################# Function for NegBin model #######################

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


compute_grad_log_full_negbin_geom <- function(s_hat, alpha_bar_hat, nstar, p,
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


compute_log_ratio_full_negbin_geom <- function(s_hat_prop, alpha_bar_hat_prop,
                                               s_hat_curr, alpha_bar_hat_curr, 
                                               nstar, p, n, K, counts,  
                                               a_s, b_s, a_alpha, b_alpha){
  
  es_prop <- exp(s_hat_prop)
  ea_bar_prop <- exp(alpha_bar_hat_prop)
  es_curr <- exp(s_hat_curr)
  ea_bar_curr <- exp(alpha_bar_hat_curr)
  
  p1_prop <- exp(lgamma(n+es_prop) + lgamma(es_prop + ea_bar_prop) -
                   lgamma(es_prop) - lgamma(es_prop + ea_bar_prop +n) )
  p1_curr <- exp(lgamma(n+es_curr) + lgamma(es_curr + ea_bar_curr) -
                   lgamma(es_curr) - lgamma(es_curr + ea_bar_curr +n) )
  
  v <- lgamma(ea_bar_prop + counts) - lgamma(1 + ea_bar_prop) + lgamma(es_prop+n - counts) -
    lgamma(es_prop) - lgamma(ea_bar_curr + counts) + lgamma(1 + ea_bar_curr) - 
    lgamma(es_curr+n - counts) + lgamma(es_curr)
  
  res <- K*(lgamma(n+es_curr+ea_bar_curr) + lgamma(es_prop + ea_bar_prop) -
              lgamma(es_curr + ea_bar_curr) - lgamma(es_prop+ea_bar_prop+n)) -
    (nstar + K)*(log(1-(1-p)*p1_prop) - log(1-(1-p)*p1_curr)) + 
    sum(v) + 
    (K+a_alpha)*(alpha_bar_hat_prop - alpha_bar_hat_curr) - b_alpha*(ea_bar_prop - ea_bar_curr) + 
    a_s*(s_hat_prop - s_hat_curr) - b_s*(es_prop - es_curr) 
  
  
  return(res)
  
}


compute_log_ratio_q_negbin_geom <- function(s_hat_prop, alpha_bar_hat_prop,
                                            s_hat_curr, alpha_bar_hat_curr, 
                                            tau, grad_log_full_curr, grad_log_full_prop){
  
  norm_prop <- norm(c(s_hat_curr, alpha_bar_hat_curr) - 
                      c(s_hat_prop, alpha_bar_hat_prop) -
                      tau*grad_log_full_prop, type="2")^2
  
  norm_curr <- norm( c(s_hat_prop, alpha_bar_hat_prop) - 
                       c(s_hat_curr, alpha_bar_hat_curr) -
                       tau*grad_log_full_curr, type="2")^2
  
  res <- 1/(4*tau)*(norm_curr - norm_prop)
  
  return (res)
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

gibbs_sampler_negbin_geometric_prior_pars <- function(Z, 
                                                      nstar_0, p_0, s_0, alpha_bar_0,
                                                      q_star, alpha_p, beta_p, a_alpha, b_alpha, a_s, b_s,
                                                      tau, fixed,
                                                      S, n_burnin, thin, seed){
  
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
  nstar_vec <- vector(length = number_saved_iterations)
  p_vec <- vector(length = number_saved_iterations)
  alpha_bar_vec <- vector(length = number_saved_iterations)
  s_vec <- vector(length = number_saved_iterations)
  
  # Set initial values
  nstar <- nstar_0
  p <- p_0
  alpha_bar <- alpha_bar_0
  s <- s_0
  
  # index saved iterations (after burn-in and thinning satisfied)
  l <- 1
  
  for (q in 1:S){
    
    t <- lgamma(s + n) + lgamma(s + alpha_bar) - lgamma(s) - lgamma(s+alpha_bar+n)
    
    ################################################################
    ############# Draw p | Z, alpha_bar, s, nstar ###################
    ###############################################################
    
    if (fixed[2] == FALSE){
      # Update Ntilde | p, Z, alpha_bar, s, nstar
      
      Ntilde <- K + rnbinom(1, size = nstar + K, prob = 1 - (1-p)* exp(t))
      
      # Update p | Ntilde, Z, alpha_bar, s, nstar
      
      p <- rbeta(1, shape1 = nstar + alpha_p, shape2 = Ntilde + beta_p)
      
    }
    
    ################################################################
    ############# Draw nstar | Z, alpha_bar, s, p ###################
    ###############################################################
    
    if (fixed[1] == FALSE){
      nstar <- 1 + rnbinom(1, size = K +1, 1 - ((1-q_star)*p)/(1 - (1-p)*exp(t)) )
    }
    
    ################################################################
    ############# Draw ( s, alpha_bar) | Z, p, nstar ##################
    ###############################################################
    
    # In order to update s, alpha_bar, we update s_hat, alpha_bar_hat, 
    # defined as the logarithm of s, alpha_bar
    
    ### Current values for s_hat, alpha_bar_hat
    s_hat_curr <- log(s)
    alpha_bar_hat_curr <- log(alpha_bar)
    
    ### Propose values for s_hat, alpha_bar_hat
    # Compute the gradient of log-full conditional for MALA 
    # (the order of returned variables is: s_hat, alpha_bar_hat)
    grad_log_full_curr <- compute_grad_log_full_negbin_geom(s_hat_curr, alpha_bar_hat_curr, 
                                                            nstar, p, n, K, counts, 
                                                            a_s, b_s, a_alpha, b_alpha)
    
    if (fixed[3] == TRUE){
      s_hat_prop <- s_hat_curr
    } else {
      s_hat_prop <- s_hat_curr + tau*grad_log_full_curr[1] + sqrt(2*tau)*rnorm(1)
    }
    
    if (fixed[4] == TRUE){
      alpha_bar_hat_prop <- alpha_bar_hat_curr
    } else {
      alpha_bar_hat_prop <- alpha_bar_hat_curr + tau*grad_log_full_curr[2] + sqrt(2*tau)*rnorm(1)
    }
    
    ### Acceptance probability 
    # Compute the log ratio of the full-cond in prop point and curr point
    log_ratio_full <- compute_log_ratio_full_negbin_geom(s_hat_prop, alpha_bar_hat_prop,
                                                         s_hat_curr, alpha_bar_hat_curr, 
                                                         nstar, p, n, K, counts,   
                                                         a_s, b_s, a_alpha, b_alpha)
    
    # Compute the log ratio of the terms related to the proposal q
    grad_log_full_prop <- compute_grad_log_full_negbin_geom(s_hat_prop, alpha_bar_hat_prop, 
                                                            nstar, p, n, K, counts, 
                                                            a_s, b_s, a_alpha, b_alpha)
    
    log_ratio_q <- compute_log_ratio_q_negbin_geom(s_hat_prop, alpha_bar_hat_prop,
                                                   s_hat_curr, alpha_bar_hat_curr, tau,
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
      
      nstar_vec[l] <- nstar
      p_vec[l] <- p
      alpha_bar_vec[l] <- alpha_bar
      s_vec[l] <- s
      
      l <- l+1
    }
    
  }
  
  nstar_vec <- nstar_vec[1:(l-1)]
  p_vec <- p_vec[1:(l-1)]
  alpha_bar_vec <- alpha_bar_vec[1:(l-1)]
  s_vec <- s_vec[1:(l-1)]
  
  return (list("nstar_vec" = nstar_vec, "p_vec" = p_vec, 
               "alpha_bar_vec" = alpha_bar_vec, "s_vec" = s_vec))
  
}

gibbs_sampler_negbin_geometric_fixed_pars <- function(Z, 
                                                      s_0, alpha_bar_0,
                                                      nstar, p, a_alpha, b_alpha, a_s, b_s,
                                                      tau, 
                                                      S, n_burnin, thin, seed){
  
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
    grad_log_full_curr <- compute_grad_log_full_negbin_geom(s_hat_curr, alpha_bar_hat_curr, 
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
    grad_log_full_prop <- compute_grad_log_full_negbin_geom(s_hat_prop, alpha_bar_hat_prop, 
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
  
  return (list("alpha_bar_vec" = alpha_bar_vec, "s_vec" = s_vec))
  
}



generate_Kmn_chain_negbin <- function(nstar, p, alpha_chain, theta_chain, M, n, Kn = 0){
  
  if (n == 0 & Kn != 0){
    stop("if n=0, the number of observed features Kn must be 0!")
  }
  
  S <- length(alpha_chain)
  
  if (length(nstar) == 1){
    nstar_chain <- rep(nstar, S)
  } else {
    nstar_chain <- nstar
  }
  
  if (length(p) == 1){
    p_chain <- rep(p, S)
  } else {
    p_chain <- p
  }
  
  M_vec <- 1:M
  kmn_chain <- matrix(NA, nrow = M, ncol = S )
  for (q in 1:S){
    nstar <- nstar_chain[q]
    p <- p_chain[q]
    alpha <- alpha_chain[q]
    theta <- theta_chain[q]
    
    par_0 = (1-p)*exp(lgamma(theta) - lgamma(theta+alpha))
    par_1 = exp(lgamma( theta + alpha +n)- lgamma(theta +n))
    par_2 = exp(lgamma(theta+alpha+n+M_vec) - lgamma(theta+n+M_vec))
    
    p_bar <- 1- par_0 *(par_1 - par_2)/(1-par_0*par_2)
    
    kmn_chain[,q] <- rnbinom(M, nstar + Kn, p_bar)
  }
  
  return (kmn_chain)
}



generate_Ntilde_chain_negbin <- function(nstar, p,
                                         alpha_chain_negbin, theta_chain_negbin,
                                         n, Kn){
  
  S <- length(alpha_chain_negbin)
  
  if (length(nstar) == 1){
    nstar_chain_negbin <- rep(nstar, S)
  } else {
    nstar_chain_negbin <- nstar
  }
  
  if (length(p) == 1){
    p_chain_negbin <- rep(p, S)
  } else {
    p_chain_negbin <- p
  }  
  
  Ntilde_chain <- vector(length = S )
  for (q in 1:S){
    nstar <- nstar_chain_negbin[q]
    p <- p_chain_negbin[q]
    alpha <- alpha_chain_negbin[q]
    theta <- theta_chain_negbin[q]
    
    negbin_par <- 1 - (1 - p)*exp(lgamma(theta+alpha+n) - lgamma(theta+alpha) - 
                                    lgamma(theta+n) + lgamma(theta))
    
    Ntilde_chain[q] <- Kn + rnbinom(1, nstar+Kn, negbin_par)
  }
  
  return (Ntilde_chain)
}


####################### Functions for Gamma IBP model #############

cppFunction('std::vector<double> stable_sum_M_all_gamma_IBP(double alpha,double theta,int m, int n){
  
  std::vector<double> sum_M; 
  sum_M.resize(m);
  
  sum_M[0] = exp(lgamma(alpha + theta + n) -
    lgamma(alpha + theta) -
    lgamma(theta + n + 1) +
    lgamma(theta + 1) )  ;
    
  for (int j=2; j < m+1; j++){
    double par = 0;
    for (int h =1; h<j+1; h++){
      par += exp(lgamma(alpha + theta + n + h - 1) -
        lgamma(alpha + theta) -
        lgamma(theta + n + h) +
        lgamma(theta + 1) ) ;
    }
    sum_M[j-1] = par;
  }
  
  return  sum_M;
   
}')

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

compute_grad_log_full_gamma_ibp <- function(s_hat, alpha_hat, Gam,
                                            n, K, counts, 
                                            a_s, b_s, a_alpha, b_alpha){
  
  es <- exp(s_hat)
  ealpha <- exp(alpha_hat)
  
  # derivative wrt s_hat
  
  p2 <- exp( lgamma(es + (1:n) -1) - lgamma(es) - 
               lgamma(es - ealpha/(1+ealpha) + (1:n)) + lgamma(es - ealpha/(1+ealpha) + 1) )
  
  p3 <- digamma(es + (1:n) -1) - digamma(es) - 
    digamma(es - ealpha/(1+ealpha) + (1:n)) + digamma(es - ealpha/(1+ealpha) + 1)
  
  ds_hat <- a_s + es*(K*( digamma(es - ealpha/(1+ealpha) + 1 ) - digamma(es - ealpha/(1+ealpha) + n) ) -
                        Gam*sum(p2 * p3 ) +
                        sum(digamma(n - counts + es)) - K*digamma(es) - b_s)
  
  # derivative wrt alpha_hat
  dalpha_hat <- a_alpha + ealpha/((1 + ealpha)^2) *
    ( K*(digamma(es - ealpha/(1+ealpha) + n ) - digamma(es - ealpha/(1+ealpha) + 1) ) -
        Gam*sum( p2 * (digamma(es - ealpha/(1+ealpha) + (1:n)) - digamma(es - ealpha/(1+ealpha) + 1))) -
        sum(digamma(1/(1 + ealpha) + counts - 1)) + 
        K*digamma(1/(1+ealpha)) - 
        (a_alpha + b_alpha)*(1+ ealpha))
  
  return (c(ds_hat, dalpha_hat))
}


compute_log_ratio_full_gamma_ibp <- function(s_hat_prop, alpha_hat_prop,
                                             s_hat_curr, alpha_hat_curr, 
                                             Gam, n, K, counts,  
                                             a_s, b_s, a_alpha, b_alpha){
  
  es_prop <- exp(s_hat_prop)
  ealpha_prop <- exp(alpha_hat_prop)
  es_curr <- exp(s_hat_curr)
  ealpha_curr <- exp(alpha_hat_curr)
  
  p2_prop <- exp( lgamma(es_prop + (1:n) -1) - lgamma(es_prop) - 
                    lgamma(es_prop - ealpha_prop/(1+ealpha_prop) + (1:n)) + 
                    lgamma(es_prop - ealpha_prop/(1+ealpha_prop) + 1) )
  
  p2_curr <- exp( lgamma(es_curr + (1:n) -1) - lgamma(es_curr) - 
                    lgamma(es_curr - ealpha_curr/(1+ealpha_curr) + (1:n)) + 
                    lgamma(es_curr - ealpha_curr/(1+ealpha_curr) + 1) )
  
  v <- lgamma(1/(1 + ealpha_prop) + counts - 1) - lgamma(1/(1 + ealpha_prop)) +
    lgamma(es_prop+ n - counts) - lgamma(es_prop) -
    lgamma(1/(1 + ealpha_curr) + counts - 1) + lgamma(1/(1 + ealpha_curr)) - 
    lgamma(es_curr+n - counts) + lgamma(es_curr)
  
  res <- K*(lgamma(es_prop - ealpha_prop/(1+ealpha_prop) + 1) - 
              lgamma(es_prop - ealpha_prop/(1+ealpha_prop) + n) -
              lgamma(es_curr - ealpha_curr/(1+ealpha_curr) + 1) + 
              lgamma(es_curr - ealpha_curr/(1+ealpha_curr) + n)) -
    Gam*sum(p2_prop - p2_curr) +
    sum(v) + 
    a_alpha*(alpha_hat_prop - alpha_hat_curr) - 
    (a_alpha + b_alpha)*(log(1 + ealpha_prop) - log(1 + ealpha_curr)) + 
    a_s*(s_hat_prop - s_hat_curr) -
    b_s*(es_prop - es_curr)
  
  return(res)
  
}


compute_log_ratio_q_gamma_ibp <- function(s_hat_prop, alpha_hat_prop,
                                          s_hat_curr, alpha_hat_curr, 
                                          tau, grad_log_full_curr, grad_log_full_prop){
  
  norm_prop <- norm(c(s_hat_curr, alpha_hat_curr) - 
                      c(s_hat_prop, alpha_hat_prop) -
                      tau*grad_log_full_prop, type="2")^2
  
  norm_curr <- norm( c(s_hat_prop, alpha_hat_prop) - 
                       c(s_hat_curr, alpha_hat_curr) -
                       tau*grad_log_full_curr, type="2")^2
  
  res <- 1/(4*tau)*(norm_curr - norm_prop)
  
  return (res)
}


gibbs_sampler_gamma_ibp <- function(Z, 
                                    a_0, b_0, s_0, alpha_0,
                                    p, r, t, a_alpha, b_alpha, a_s, b_s,
                                    sigq_s, sigq_alpha, fixed,
                                    S, n_burnin, thin, seed){
  
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
    
    if (fixed[1] == FALSE){
      
      # Update a | Gam, Z, alpha, s, b
      
      a <- 1 + rpois(1, b*Gam*(1-p) )
      
    }
    
    if (fixed[2] == FALSE){
      
      # Update b | Gam, Z, alpha, s, a
      
      b <- rgamma(1, shape = a + r, rate = Gam + t )
      
    }
    
    
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
    if (fixed[4] == TRUE){
      s_hat_prop <- s_hat_curr
    } else {
      s_hat_prop <- rnorm(n=1, mean = s_hat_curr, sd = sqrt(sigq_s) ) 
    }
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
    if (fixed[3] == TRUE){
      alpha_hat_prop <- alpha_hat_curr
    } else {
      alpha_hat_prop <- rnorm(n=1, mean = alpha_hat_curr, sd = sqrt(sigq_alpha) )
    }
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
  
  return (list("a_vec" = a_vec, "b_vec" = b_vec, 
               "alpha_vec" = alpha_vec, "s_vec" = s_vec))
  
}


generate_Kmn_chain_gamma_ibp <- function(a_chain, b_chain, alpha_chain, theta_chain, M, n, Kn = 0){
  
  if (n == 0 & Kn != 0){
    stop("if n=0, the number of observed features Kn must be 0!")
  }
  
  S <- length(a_chain)
  M_vec <- 1:M
  kmn_chain <- matrix(NA, nrow = M, ncol = S )
  for (q in 1:S){
    a <- a_chain[q]
    b <- b_chain[q]
    alpha <- alpha_chain[q]
    theta <- theta_chain[q]
    
    if (n == 0){
      gamma_a_t_n <- 0
    } else {
      gamma_a_t_n <- sum(exp(lgamma(alpha + theta + (1:n) - 1) - lgamma(alpha + theta) -
                               lgamma(theta + (1:n)) + lgamma(theta +1) ) )
    }
    
    # sum_M <- sapply(M_vec, function(m) sum(exp(lgamma(alpha + theta + n + (1:m) - 1) -
    #                                              lgamma(alpha + theta) -
    #                                              lgamma(theta + n + (1:m)) +
    #                                              lgamma(theta +1) ) ) )
    sum_M <- stable_sum_M_all_gamma_IBP(alpha, theta, M, n)
    #print(sum_M)
    kmn_chain[,q] <- rnbinom(M, a + Kn, (gamma_a_t_n + b)/(gamma_a_t_n + b + sum_M))
    #print(paste0("iteration: ", q))
    # pbar <- p_kmn_all_gamma_IBP(alpha, theta, M, n, b)
    # print(paste0("pbar: ", pbar))
    # kmn_chain[,q] <- rnbinom(M, a + Kn, pbar)
  }
  
  return (kmn_chain)
}


#################### Functions for SB-SP model #################
compute_der_log_full_sb_sp <- function(alpha_hat, 
                                       c, beta, Gam, n, K, counts,
                                       a_alpha, b_alpha){
  
  ealpha <- exp(alpha_hat)
  
  p2 <- exp( lgamma((1:n)) - 
               lgamma(1 - ealpha/(1+ealpha) + (1:n)) + lgamma(2 - ealpha/(1+ealpha) ) )
  
  # contribution of the Gamma-IBP
  dalpha_hat <- a_alpha + ealpha/((1 + ealpha)^2) *
    ( K*(digamma(1 - ealpha/(1+ealpha) + n ) - digamma(2 - ealpha/(1+ealpha) ) ) -
        Gam*sum( p2 * (digamma(1 - ealpha/(1+ealpha) + (1:n)) - digamma(2 - ealpha/(1+ealpha) ))) -
        sum(digamma(1/(1 + ealpha) + counts - 1)) + 
        K*digamma(1/(1+ealpha)) - 
        (a_alpha + b_alpha)*(1+ ealpha))
  
  # due to the fact alpha is in the prior on Gam
  dalpha_hat <- dalpha_hat - (c + 1) + Gam*beta/ealpha
  
  return (dalpha_hat)
  
}


compute_log_ratio_full_sb_sp <- function(alpha_hat_prop,
                                         alpha_hat_curr, 
                                         c, beta, Gam, n, K, counts,
                                         a_alpha, b_alpha){
  
  ealpha_prop <- exp(alpha_hat_prop)
  ealpha_curr <- exp(alpha_hat_curr)
  
  p2_prop <- exp( lgamma( (1:n) ) - 
                    lgamma(1 - ealpha_prop/(1+ealpha_prop) + (1:n)) + 
                    lgamma(2 - ealpha_prop/(1+ealpha_prop) ) )
  
  p2_curr <- exp( lgamma( (1:n) ) - 
                    lgamma(1 - ealpha_curr/(1+ealpha_curr) + (1:n)) + 
                    lgamma(2 - ealpha_curr/(1+ealpha_curr) ) )
  
  v <- lgamma(1/(1 + ealpha_prop) + counts - 1) - lgamma(1/(1 + ealpha_prop)) -
    lgamma(1/(1 + ealpha_curr) + counts - 1) + lgamma(1/(1 + ealpha_curr)) 
  
  # log-ratio of the Gamma-IBP contribution
  res <- K*(lgamma(2 - ealpha_prop/(1+ealpha_prop) ) - 
              lgamma(1 - ealpha_prop/(1+ealpha_prop) + n) -
              lgamma(2 - ealpha_curr/(1+ealpha_curr) ) + 
              lgamma(1 - ealpha_curr/(1+ealpha_curr) + n)) -
    Gam*sum(p2_prop - p2_curr) +
    sum(v) + 
    a_alpha*(alpha_hat_prop - alpha_hat_curr) - 
    (a_alpha + b_alpha)*(log(1 + ealpha_prop) - log(1 + ealpha_curr))
  
  # due to the fact alpha is in prior on Gam
  res <- res + (c+1)*(alpha_hat_curr - alpha_hat_prop) + 
    Gam*beta*(1/ealpha_curr - 1/ealpha_prop)
  
  
  return (res)
  
}


compute_log_ratio_q_sb_sp <- function(alpha_hat_prop,
                                      alpha_hat_curr, tau,
                                      der_log_full_curr, der_log_full_prop){
  
  norm_prop <- (alpha_hat_curr - alpha_hat_prop - tau*der_log_full_prop)^2
  
  norm_curr <- (alpha_hat_prop - alpha_hat_curr - tau*der_log_full_curr)^2
  
  res <- 1/(4*tau)*(norm_curr - norm_prop)
  
  return (res)
  
}


gibbs_sampler_sb_sp <- function(Z, 
                                c_0, beta_0, alpha_0,
                                p, r, t, a_alpha, b_alpha, 
                                tau, fixed,
                                S, n_burnin, thin, seed){
  
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
  c_vec <- vector(length = number_saved_iterations)
  beta_vec <- vector(length = number_saved_iterations)
  alpha_vec <- vector(length = number_saved_iterations)
  
  # Set initial values
  c <- c_0
  beta <- beta_0
  alpha <- alpha_0
  
  # index saved iterations (after burn-in and thinning satisfied)
  l <- 1
  
  for (q in 1:S){
    
    ################################################################
    ################# Draw c, beta, alpha | Z   ###################
    ###############################################################
    
    # Update Gam | Z, alpha, c, beta
    
    gamma_a_1ma_n <- sum(exp(lgamma((1:n)) - 
                               lgamma(1 - alpha + (1:n)) + lgamma(2 - alpha) ) )
    
    Gam <- rgamma(1, shape = K + c + 1, rate = beta*(1-alpha)/alpha + gamma_a_1ma_n)
    
    if (fixed[1] == FALSE){
      
      # Update c | Gam, Z, alpha, beta
      
      c <- rpois(1, beta*Gam*(1-p)*(1-alpha)/alpha )
      
    }
    
    if (fixed[2] == FALSE){
      
      # Update beta | Gam, Z, alpha, c
      
      beta <- rgamma(1, shape = c + 1 + r, rate = Gam*(1-alpha)/alpha + t )
      
    }
    
    
    ################################################################
    ############# Draw alpha | Z, Gam, c, beta ##################
    ###############################################################
    
    # In order to update alpha, we update alpha_hat, 
    # defined as transformation of alpha
    
    ### Current values for alpha_hat
    alpha_hat_curr <- log(alpha/(1-alpha))
    
    ### Propose values for alpha_hat
    # Compute the gradient of log-full conditional for MALA 
    der_log_full_curr <- compute_der_log_full_sb_sp(alpha_hat_curr, 
                                                    c, beta, Gam, n, K, counts,
                                                    a_alpha, b_alpha)
    
    if (fixed[3] == TRUE){
      alpha_hat_prop <- alpha_hat_curr
    } else {
      alpha_hat_prop <- alpha_hat_curr + tau*der_log_full_curr + sqrt(2*tau)*rnorm(1)
    }
    
    ### Acceptance probability 
    # Compute the log ratio of the full-cond in prop point and curr point
    log_ratio_full <- compute_log_ratio_full_sb_sp(alpha_hat_prop,
                                                   alpha_hat_curr, 
                                                   c, beta, Gam, n, K, counts,
                                                   a_alpha, b_alpha)
    
    # Compute the log ratio of the terms related to the proposal q
    der_log_full_prop <- compute_der_log_full_sb_sp(alpha_hat_prop,
                                                    c, beta, Gam, n, K, counts,
                                                    a_alpha, b_alpha)
    
    log_ratio_q <- compute_log_ratio_q_sb_sp(alpha_hat_prop,
                                             alpha_hat_curr, tau,
                                             der_log_full_curr, der_log_full_prop)
    
    # Compute acceptance probability
    log_acc_prob <- log_ratio_full + log_ratio_q
    acc_prob <- min(1, exp(log_acc_prob))
    
    # Decide if accept or not the new alpha
    if (runif(1) < acc_prob){ # accept
      alpha <- exp(alpha_hat_prop)/(1 + exp(alpha_hat_prop))
    }
    
    
    #######################################################################
    
    # Store parameters if burn-in is over and once every "thin" iteration
    if ((q > n_burnin) & (q %% thin == 0) ){
      print(paste0("iteration: ", q))
      
      c_vec[l] <- c
      beta_vec[l] <- beta
      alpha_vec[l] <- alpha
      
      l <- l+1
    }
    
  }
  
  c_vec <- c_vec[1:(l-1)]
  beta_vec <- beta_vec[1:(l-1)]
  alpha_vec <- alpha_vec[1:(l-1)]
  
  return (list("c_vec" = c_vec, "beta_vec" = beta_vec, 
               "alpha_vec" = alpha_vec))
  
}



####
##### MODEL 4: the log-normal model ############################
####

library(tidyverse)

# set number of individuals
Ns <- c(20, 40, 80) 
L <- 300

# maximum number of features
H <- 500

# vector of ai's 
seed = 123456
set.seed(seed)
as <- exp(rnorm(H, 0, 1))

# choose c such that c*as <= 0.5
c <- 1 / max(as)
# define pi's
pis <- c*as

###### 1) Set parameters for the 3 models ###############

# Set desired value of E[N] = Nbar
Nbars <- c(200, 400, 600)
# Set value of c_fr such that Var(N)= c_fr * Nbar, when two or more parameters
c_fr <- 10

########### 1.1) Set parameters for BB (Poisson and Neg-Bin) for alpha, theta

# Set hyperparameters
a_alpha_bb <- 1
b_alpha_bb <- 0.1
print(paste0("E(alpha_bar) = ", a_alpha_bb/ b_alpha_bb))
print(paste0("Var(alpha_bar) = ", a_alpha_bb/ (b_alpha_bb^2)))
a_s_bb <- 2
b_s_bb <- 0.2
print(paste0("E(s) = ", a_s_bb/ b_s_bb))
print(paste0("Var(s) = ", a_s_bb/ (b_s_bb^2)))

# # Set initial values for NB (prior on parameters)
# nstar_0_nb <- 100
# p_0_nb <- 0.2
# Set initial values for other parameters
alpha_bar_0_bb <- 1
s_0_bb <- 1


########## 1.2) Set parameters for the Gamma IBP and SB-SP

# Set prior hyperparameters for Gamma IBP
p_ibp <- 0.05
print(paste0("E(a) = ", 1/ p_ibp))
print(paste0("Var(a) = ", (1-p_ibp)/ (p_ibp^2)) )
r_ibp <- 1
t_ibp <- 0.1
print(paste0("E(b) = ", r_ibp/ t_ibp))
print(paste0("Var(b) = ", r_ibp/ (t_ibp^2)))
a_alpha_ibp <- 2
b_alpha_ibp <- 2
print(paste0("E(alpha) = ", a_alpha_ibp/ (a_alpha_ibp + b_alpha_ibp)))
print(paste0("Var(alpha) = ", a_alpha_ibp*b_alpha_ibp/ (a_alpha_ibp + b_alpha_ibp)^2 /(a_alpha_ibp+b_alpha_ibp+1)))
a_s_ibp <- 2
b_s_ibp <- 0.2
print(paste0("E(s) = ", a_s_ibp/ b_s_ibp))
print(paste0("Var(s) = ", a_s_ibp/ (b_s_ibp^2)))

print(paste0("E(theta) = ", a_s_ibp/ b_s_ibp -  a_alpha_ibp/ (a_alpha_ibp + b_alpha_ibp)))
print(paste0("Var(theta) = ", a_s_ibp/ (b_s_ibp^2) + 
               a_alpha_ibp*b_alpha_ibp/ (a_alpha_ibp + b_alpha_ibp)^2 /(a_alpha_ibp+b_alpha_ibp+1)))


# # Set prior hyperparameters for SB-SP
# p_sp <- 0.08
# print(paste0("E(c) = ", (1-p_sp)/ p_sp))
# print(paste0("Var(c) = ", (1-p_sp)/ (p_sp^2)) )
# r_sp <- 0.1
# t_sp <- 0.01
# print(paste0("E(beta) = ", r_sp/ t_sp))
# print(paste0("Var(beta) = ", r_sp/ (t_sp^2)))
# a_alpha_sp <- 2
# b_alpha_sp <- 2
# print(paste0("E(alpha) = ", a_alpha_sp/ (a_alpha_sp + b_alpha_sp)))
# print(paste0("Var(alpha) = ", a_alpha_sp*b_alpha_sp/ (a_alpha_sp + b_alpha_sp)^2 /(a_alpha_sp+b_alpha_sp+1)))


# Set initial values for the parameters of Gamma IBP
a_0_ibp <- 5
b_0_ibp <- 1
alpha_0_ibp <- 0.5
s_0_ibp <- 15

# # Set initial values for the parameters of SB-SP
# c_0_sp <- 10
# beta_0_sp <- 10
# alpha_0_sp <- 0.5



##### Single dataset -> Ntilde (Poiss/NB) and extrapolation (Poiss/NB/Gamma) #####

set.seed(seed)
data_mat <- matrix(rbinom(L*H, size = 1, prob = rep(pis, L)),
                   nrow = L, ncol = H, byrow = T )
data_list <- create_features_list(data_mat)
#plot_trajectory(data_list)

num_feat <- vector(length = length(Ns))
for (j in 1:length(Ns)){
  N <- Ns[j]
  num_feat[j] <- sum(colSums(data_mat[1:N, ]) > 0)
}
#print("N. of observed features in the sample: ")
#print(num_feat)

########## Set MCMC parameters (common to all 3 models)

S_poiss <- S_negbin <- 3*10^4
S_ibp <-  5*10^4
n_burnin_poiss <- n_burnin_negbin <- n_burnin_ibp<-  5*10^3
thin_poiss <- thin_negbin <- thin_ibp <- 2
seed <- 1234
number_saved_iterations_poiss <- (S_poiss - n_burnin_poiss)/thin_poiss
number_saved_iterations_negbin <- (S_negbin - n_burnin_negbin)/thin_negbin
number_saved_iterations_ibp <- (S_ibp - n_burnin_ibp)/thin_ibp

###### 2) Run the algorithms ###############
labels_comb <- paste(rep(paste("N", Ns, sep = "."), each = length(Nbars)+1),
                     c(paste("Nbar", Nbars, sep = "."),"Nbar.emp"), sep=":")

gg_ntilde_poiss <- data.frame(matrix(nrow = number_saved_iterations_poiss, ncol = length(Ns)*(length(Nbars)+1)))
colnames(gg_ntilde_poiss) <- labels_comb
gg_ntilde_negbin <- data.frame(matrix(nrow = number_saved_iterations_negbin, ncol = length(Ns)*(length(Nbars)+1)))
colnames(gg_ntilde_negbin) <- labels_comb
# gg_ntilde_negbin_prior <- data.frame(matrix(nrow = number_saved_iterations_negbin, ncol = length(Ns)))
# colnames(gg_ntilde_negbin_prior) <- paste("N", Ns, sep = ".")

list_kmn_pred_test_poiss <- vector(mode="list", length = length(Ns)*(length(Nbars)+1))
names(list_kmn_pred_test_poiss) <- labels_comb
list_kmn_pred_test_negbin <- vector(mode="list", length = length(Ns)*(length(Nbars)+1))
names(list_kmn_pred_test_negbin) <- labels_comb
# list_kmn_pred_test_negbin_prior <- vector(mode="list", length = length(Ns))
# names(list_kmn_pred_test_negbin_prior) <- paste("N", Ns, sep = ".")
list_kmn_pred_test_ibp <- vector(mode="list", length = length(Ns)*(length(Nbars)+1))
names(list_kmn_pred_test_ibp) <- labels_comb
# list_kmn_pred_test_sp <- vector(mode="list", length = length(Ns))
# names(list_kmn_pred_test_sp) <- paste("N", Ns, sep = ".")

params_poiss <- data.frame(matrix(nrow = number_saved_iterations_poiss, ncol = 2*length(Ns)*(length(Nbars)+1)))
colnames(params_poiss) <- paste(c("alpha", "theta"), rep(labels_comb, each = 2), sep = ":")
params_negbin <- data.frame(matrix(nrow = number_saved_iterations_negbin, ncol = 2*length(Ns)*(length(Nbars)+1)))
colnames(params_negbin) <- paste(c("alpha", "theta"), rep(labels_comb, each = 2), sep = ":")
# params_negbin_prior <- data.frame(matrix(nrow = number_saved_iterations_negbin, ncol = 4*length(Ns)))
# colnames(params_negbin_prior) <- paste(c("nstar", "p", "alpha", "theta"), rep(Ns, each = 4), sep = ".")
params_ibp <- data.frame(matrix(nrow = number_saved_iterations_ibp, ncol = 4*length(Ns)*(length(Nbars)+1)))
colnames(params_ibp) <- paste(c("a", "b", "alpha", "theta"), rep(labels_comb, each = 4), sep = ":")
# params_sp <- data.frame(matrix(nrow = number_saved_iterations_ibp, ncol = 3*length(Ns)))
# colnames(params_sp) <- paste(c("c", "beta", "alpha"), rep(Ns, each = 3), sep = ".")

for (j in 1:length(Ns)){
  N <- Ns[j]
  M <- L - N
  
  train_mat <- data_mat[1:N,]
  test_mat <- data_mat[(N+1):L, ]
  # convert the binary matrix into list of features
  train_list <- create_features_list(train_mat)
  test_list <- create_features_list(test_mat)
  
  Kn = ncol(train_mat[,colSums(train_mat) > 0])
  print(paste0("Number of observed features: ", Kn))
  
  ##### 2.2) Run the empirical Nbar #####
  
  Nbar_emp <- beta_binomial_estimator(train_mat)
  
  # Set prior hyperparameters specific for poisson, with fixed lambda
  lambda_poiss <- Nbar_emp
  print("Poisson params")
  print(paste0("E(N) = ", lambda_poiss))
  print(paste0("Var(N) = ", lambda_poiss))
  
  # Set prior hyperparameters specific for NB, with fixed parameters
  nstar_nb <- Nbar_emp/(c_fr - 1)
  p_nb <- 1/c_fr
  print("NB params (fixed)")
  print(paste0("E(N) = ", nstar_nb*(1-p_nb)/p_nb ))
  print(paste0("Var(N) = ", nstar_nb*(1-p_nb)/(p_nb^2) ))
  
  # Label for accessing element of structures related to N and Nbar
  lab_comb <- paste0("N.",N,":Nbar.emp")
  
  lab_alpha <- paste0("alpha:",lab_comb)
  lab_theta <- paste0("theta:",lab_comb)
  lab_a <- paste0("a:",lab_comb)
  lab_b <- paste0("b:",lab_comb)
  
  ################# 3) Run BB + Poisson ###########
  
  # Set tau for MALA
  tau_poiss <- 0.1
  
  output_poiss <- gibbs_sampler_poiss_fixed_lambda(Z = train_mat,
                                                   alpha_bar_0_bb, s_0_bb,
                                                   lambda_poiss, a_alpha_bb, b_alpha_bb, a_s_bb, b_s_bb,
                                                   tau_poiss,
                                                   S_poiss, n_burnin_poiss, thin_poiss, seed)
  
  n_saved_iter_poiss <- length(output_poiss$s_vec)
  s_chain_poiss <- output_poiss$s_vec
  alpha_bar_chain_poiss <- output_poiss$alpha_bar_vec
  alpha_chain_poiss <- - alpha_bar_chain_poiss
  theta_chain_poiss <- s_chain_poiss + alpha_bar_chain_poiss
  
  params_poiss[[lab_alpha]] <- alpha_chain_poiss
  params_poiss[[lab_theta]] <- theta_chain_poiss
  
  
  ####### 4) Run BB + Negative-Binomial (fixed parameters) ############
  
  # Set tau for MALA
  tau_nb <- 0.1
  
  output_negbin <- gibbs_sampler_negbin_geometric_fixed_pars(Z = train_mat,
                                                             s_0_bb, alpha_bar_0_bb,
                                                             nstar_nb, p_nb,
                                                             a_alpha_bb, b_alpha_bb, a_s_bb, b_s_bb,
                                                             tau_nb, 
                                                             S_negbin, n_burnin_negbin, thin_negbin, seed)
  
  n_saved_iter_negbin <- length(output_negbin$s_vec)
  s_chain_negbin <- output_negbin$s_vec
  alpha_bar_chain_negbin <- output_negbin$alpha_bar_vec
  alpha_chain_negbin <- - alpha_bar_chain_negbin
  theta_chain_negbin <- s_chain_negbin + alpha_bar_chain_negbin
  
  params_negbin[[lab_alpha]] <- alpha_chain_negbin
  params_negbin[[lab_theta]] <- theta_chain_negbin
  
  # ####### 4.b) Run BB + Negative-Binomial (prior on parameters) #######
  # 
  # # Set tau for MALA
  # tau_nb_prior <- 0.001
  # 
  # output_negbin_prior <- gibbs_sampler_negbin_geometric_prior_pars(Z = train_mat,
  #                                                            nstar_0_nb, p_0_nb, s_0_bb, alpha_bar_0_bb,
  #                                                            q_star_nb, alpha_p_nb, beta_p_nb,
  #                                                            a_alpha_bb, b_alpha_bb, a_s_bb, b_s_bb,
  #                                                            tau_nb_prior, fixed = c(F,F,F,F),
  #                                                            S_negbin, n_burnin_negbin, thin_negbin, seed)
  # 
  # n_saved_iter_negbin_prior <- length(output_negbin_prior$nstar_vec)
  # nstar_chain_negbin_prior <- output_negbin_prior$nstar_vec
  # p_chain_negbin_prior <- output_negbin_prior$p_vec
  # s_chain_negbin_prior <- output_negbin_prior$s_vec
  # alpha_bar_chain_negbin_prior <- output_negbin_prior$alpha_bar_vec
  # alpha_chain_negbin_prior <- - alpha_bar_chain_negbin_prior
  # theta_chain_negbin_prior <- s_chain_negbin_prior + alpha_bar_chain_negbin_prior
  # 
  # params_negbin_prior[[paste0("nstar.",N)]] <- nstar_chain_negbin_prior
  # params_negbin_prior[[paste0("p.",N)]] <- p_chain_negbin_prior
  # params_negbin_prior[[paste0("alpha.",N)]] <- alpha_chain_negbin_prior
  # params_negbin_prior[[paste0("theta.",N)]] <- theta_chain_negbin_prior
  # 
  ############ 5) Run IBP + Gamma #################
  
  # Set tau for MALA
  sigq_s <- 0.1
  sigq_alpha <- 0.1
  
  output_ibp <- gibbs_sampler_gamma_ibp(Z = train_mat,
                                        a_0_ibp, b_0_ibp, s_0_ibp, alpha_0_ibp,
                                        p_ibp, r_ibp, t_ibp, a_alpha_ibp, b_alpha_ibp, a_s_ibp, b_s_ibp,
                                        sigq_s, sigq_alpha, fixed = c(F,F,F,F),
                                        S_ibp, n_burnin_ibp, thin_ibp, seed)
  
  n_saved_iter_ibp <- length(output_ibp$a_vec)
  a_chain_ibp <- output_ibp$a_vec
  b_chain_ibp <- output_ibp$b_vec
  s_chain_ibp <- output_ibp$s_vec
  alpha_chain_ibp <- output_ibp$alpha_vec
  theta_chain_ibp <- s_chain_ibp - alpha_chain_ibp
  
  params_ibp[[lab_a]] <- a_chain_ibp
  params_ibp[[lab_b]] <- b_chain_ibp
  params_ibp[[lab_alpha]] <- alpha_chain_ibp
  params_ibp[[lab_theta]] <- theta_chain_ibp
  
  # ############ 5.1) Run SB-SP #################
  # 
  # # Set tau for MALA
  # tau_sp <- 0.0001
  # 
  # output_sp <- gibbs_sampler_sb_sp(Z = train_mat,
  #                                  c_0_sp, beta_0_sp, alpha_0_sp,
  #                                  p_sp, r_sp, t_sp, a_alpha_sp, b_alpha_sp,
  #                                  tau_sp, fixed = c(F,F,F),
  #                                  S_ibp, n_burnin_ibp, thin_ibp, seed)
  # 
  # n_saved_iter_sp <- length(output_sp$c_vec)
  # c_chain_sp <- output_sp$c_vec
  # beta_chain_sp <- output_sp$beta_vec
  # alpha_chain_sp <- output_sp$alpha_vec
  # 
  # params_sp[[paste0("c.",N)]] <- c_chain_sp
  # params_sp[[paste0("beta.",N)]] <- beta_chain_sp
  # params_sp[[paste0("alpha.",N)]] <- alpha_chain_sp
  
  ####### 6) Estimate limit distributions (Poiss/NB) ################
  ntilde_chain_poiss <- generate_Ntilde_chain_poiss(lambda_poiss, alpha_chain_poiss,
                                                    theta_chain_poiss, n = N, Kn)
  
  ntilde_chain_negbin <- generate_Ntilde_chain_negbin(nstar_nb, p_nb,
                                                      alpha_chain_negbin,
                                                      theta_chain_negbin, n = N, Kn)
  
  # ntilde_chain_negbin_prior <- generate_Ntilde_chain_negbin(nstar_chain_negbin_prior, 
  #                                                           p_chain_negbin_prior,
  #                                                           alpha_chain_negbin,
  #                                                           theta_chain_negbin, n = N, Kn)
  
  gg_ntilde_poiss[[lab_comb]] <- ntilde_chain_poiss
  gg_ntilde_negbin[[lab_comb]] <- ntilde_chain_negbin
  #gg_ntilde_negbin_prior[[paste0("N.",N)]] <- ntilde_chain_negbin_prior
  
  ######## 7) Extrapolation in the test set (Poiss/NB/Gamma/SP) ##############
  # Poisson
  kmn_chain_poiss <- generate_Kmn_chain_poiss(lambda_poiss, alpha_chain_poiss,
                                              theta_chain_poiss, M = M, n = N)
  
  est_ci_pred_poiss <- matrix(NA, nrow = M, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:M){
    est_ci_pred_poiss[m,] <- quantile(kmn_chain_poiss[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_pred_poiss <- list("medians" = est_ci_pred_poiss[,2],
                            "lbs" = est_ci_pred_poiss[,1],
                            "ubs" = est_ci_pred_poiss[,3])
  
  list_kmn_pred_test_poiss[[lab_comb]] <- est_ci_pred_poiss
  
  # Negative Binomial (fixed params)
  kmn_chain_negbin <- generate_Kmn_chain_negbin(nstar_nb, p_nb,
                                                alpha_chain_negbin, theta_chain_negbin,
                                                M = M, n = N, Kn)
  
  est_ci_pred_negbin <- matrix(NA, nrow = M, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:M){
    est_ci_pred_negbin[m,] <- quantile(kmn_chain_negbin[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_pred_negbin <- list("medians" = est_ci_pred_negbin[,2],
                             "lbs" = est_ci_pred_negbin[,1],
                             "ubs" = est_ci_pred_negbin[,3])
  
  list_kmn_pred_test_negbin[[lab_comb]] <- est_ci_pred_negbin
  
  # # Negative Binomial (prior on params)
  # kmn_chain_negbin_prior <- generate_Kmn_chain_negbin(nstar_chain_negbin_prior,
  #                                                     p_chain_negbin_prior,
  #                                                     alpha_chain_negbin, 
  #                                                     theta_chain_negbin,
  #                                                     M = M, n = N, Kn)
  # 
  # est_ci_pred_negbin_prior <- matrix(NA, nrow = M, ncol = 3)
  # # first column = lower bound
  # # second columns = medians
  # # third columns = upper bound
  # for (m in 1:M){
  #   est_ci_pred_negbin_prior[m,] <- quantile(kmn_chain_negbin_prior[m,], probs = c(0.025,0.5,0.975))
  # }
  # est_ci_pred_negbin_prior <- list("medians" = est_ci_pred_negbin_prior[,2],
  #                            "lbs" = est_ci_pred_negbin_prior[,1],
  #                            "ubs" = est_ci_pred_negbin_prior[,3])
  # 
  # list_kmn_pred_test_negbin_prior[[paste0("N.",N)]] <- est_ci_pred_negbin_prior
  
  
  # IBP + Gamma
  kmn_chain_ibp <- generate_Kmn_chain_gamma_ibp(a_chain_ibp, b_chain_ibp, alpha_chain_ibp,
                                                theta_chain_ibp, M = M, n = N, Kn)
  
  est_ci_pred_ibp <- matrix(NA, nrow = M, ncol = 3)
  # first column = lower bound
  # second columns = medians
  # third columns = upper bound
  for (m in 1:M){
    est_ci_pred_ibp[m,] <- quantile(kmn_chain_ibp[m,], probs = c(0.025,0.5,0.975))
  }
  est_ci_pred_ibp <- list("medians" = est_ci_pred_ibp[,2],
                          "lbs" = est_ci_pred_ibp[,1],
                          "ubs" = est_ci_pred_ibp[,3])
  
  list_kmn_pred_test_ibp[[lab_comb]] <- est_ci_pred_ibp
  
  # # SB-SP
  # kmn_chain_sp <- generate_Kmn_chain_gamma_ibp(a_chain = c_chain_sp + 1,
  #                                              b_chain = beta_chain_sp*(1-alpha_chain_sp)/alpha_chain_sp, 
  #                                              alpha_chain = alpha_chain_sp,
  #                                              theta_chain = 1 - alpha_chain_sp,
  #                                              M = M, n = N, Kn)
  # 
  # est_ci_pred_sp <- matrix(NA, nrow = M, ncol = 3)
  # # first column = lower bound
  # # second columns = medians
  # # third columns = upper bound
  # for (m in 1:M){
  #   est_ci_pred_sp[m,] <- quantile(kmn_chain_sp[m,], probs = c(0.025,0.5,0.975))
  # }
  # est_ci_pred_sp <- list("medians" = est_ci_pred_sp[,2],
  #                         "lbs" = est_ci_pred_sp[,1],
  #                         "ubs" = est_ci_pred_sp[,3])
  # 
  # list_kmn_pred_test_sp[[paste0("N.",N)]] <- est_ci_pred_sp
  
  
  #### 2.3) Run the Nbars on the grid #####
  for (v in 1:length(Nbars)){
    
    Nbar <- Nbars[v]
    
    # Set prior hyperparameters specific for poisson, with fixed lambda
    lambda_poiss <- Nbar
    print("Poisson params")
    print(paste0("E(N) = ", lambda_poiss))
    print(paste0("Var(N) = ", lambda_poiss))
    
    # Set prior hyperparameters specific for NB, with fixed parameters
    nstar_nb <- Nbar/(c_fr - 1)
    p_nb <- 1/c_fr
    print("NB params (fixed)")
    print(paste0("E(N) = ", nstar_nb*(1-p_nb)/p_nb ))
    print(paste0("Var(N) = ", nstar_nb*(1-p_nb)/(p_nb^2) ))
    
    # Label for accessing element of structures related to N and Nbar
    lab_comb <- paste0("N.",N,":Nbar.",Nbar)
    
    lab_alpha <- paste0("alpha:",lab_comb)
    lab_theta <- paste0("theta:",lab_comb)
    lab_a <- paste0("a:",lab_comb)
    lab_b <- paste0("b:",lab_comb)
    
    ################# 3) Run BB + Poisson ###########
    
    # Set tau for MALA
    tau_poiss <- 0.1
    
    output_poiss <- gibbs_sampler_poiss_fixed_lambda(Z = train_mat,
                                                     alpha_bar_0_bb, s_0_bb,
                                                     lambda_poiss, a_alpha_bb, b_alpha_bb, a_s_bb, b_s_bb,
                                                     tau_poiss,
                                                     S_poiss, n_burnin_poiss, thin_poiss, seed)
    
    n_saved_iter_poiss <- length(output_poiss$s_vec)
    s_chain_poiss <- output_poiss$s_vec
    alpha_bar_chain_poiss <- output_poiss$alpha_bar_vec
    alpha_chain_poiss <- - alpha_bar_chain_poiss
    theta_chain_poiss <- s_chain_poiss + alpha_bar_chain_poiss
    
    params_poiss[[lab_alpha]] <- alpha_chain_poiss
    params_poiss[[lab_theta]] <- theta_chain_poiss
    
    
    ####### 4) Run BB + Negative-Binomial (fixed parameters) ############
    
    # Set tau for MALA
    tau_nb <- 0.1
    
    output_negbin <- gibbs_sampler_negbin_geometric_fixed_pars(Z = train_mat,
                                                               s_0_bb, alpha_bar_0_bb,
                                                               nstar_nb, p_nb,
                                                               a_alpha_bb, b_alpha_bb, a_s_bb, b_s_bb,
                                                               tau_nb, 
                                                               S_negbin, n_burnin_negbin, thin_negbin, seed)
    
    n_saved_iter_negbin <- length(output_negbin$s_vec)
    s_chain_negbin <- output_negbin$s_vec
    alpha_bar_chain_negbin <- output_negbin$alpha_bar_vec
    alpha_chain_negbin <- - alpha_bar_chain_negbin
    theta_chain_negbin <- s_chain_negbin + alpha_bar_chain_negbin
    
    params_negbin[[lab_alpha]] <- alpha_chain_negbin
    params_negbin[[lab_theta]] <- theta_chain_negbin
    
    # ####### 4.b) Run BB + Negative-Binomial (prior on parameters) #######
    # 
    # # Set tau for MALA
    # tau_nb_prior <- 0.001
    # 
    # output_negbin_prior <- gibbs_sampler_negbin_geometric_prior_pars(Z = train_mat,
    #                                                            nstar_0_nb, p_0_nb, s_0_bb, alpha_bar_0_bb,
    #                                                            q_star_nb, alpha_p_nb, beta_p_nb,
    #                                                            a_alpha_bb, b_alpha_bb, a_s_bb, b_s_bb,
    #                                                            tau_nb_prior, fixed = c(F,F,F,F),
    #                                                            S_negbin, n_burnin_negbin, thin_negbin, seed)
    # 
    # n_saved_iter_negbin_prior <- length(output_negbin_prior$nstar_vec)
    # nstar_chain_negbin_prior <- output_negbin_prior$nstar_vec
    # p_chain_negbin_prior <- output_negbin_prior$p_vec
    # s_chain_negbin_prior <- output_negbin_prior$s_vec
    # alpha_bar_chain_negbin_prior <- output_negbin_prior$alpha_bar_vec
    # alpha_chain_negbin_prior <- - alpha_bar_chain_negbin_prior
    # theta_chain_negbin_prior <- s_chain_negbin_prior + alpha_bar_chain_negbin_prior
    # 
    # params_negbin_prior[[paste0("nstar.",N)]] <- nstar_chain_negbin_prior
    # params_negbin_prior[[paste0("p.",N)]] <- p_chain_negbin_prior
    # params_negbin_prior[[paste0("alpha.",N)]] <- alpha_chain_negbin_prior
    # params_negbin_prior[[paste0("theta.",N)]] <- theta_chain_negbin_prior
    # 
    ############ 5) Run IBP + Gamma #################
    
    # Set tau for MALA
    sigq_s <- 0.1
    sigq_alpha <- 0.1
    
    output_ibp <- gibbs_sampler_gamma_ibp(Z = train_mat,
                                          a_0_ibp, b_0_ibp, s_0_ibp, alpha_0_ibp,
                                          p_ibp, r_ibp, t_ibp, a_alpha_ibp, b_alpha_ibp, a_s_ibp, b_s_ibp,
                                          sigq_s, sigq_alpha, fixed = c(F,F,F,F),
                                          S_ibp, n_burnin_ibp, thin_ibp, seed)
    
    n_saved_iter_ibp <- length(output_ibp$a_vec)
    a_chain_ibp <- output_ibp$a_vec
    b_chain_ibp <- output_ibp$b_vec
    s_chain_ibp <- output_ibp$s_vec
    alpha_chain_ibp <- output_ibp$alpha_vec
    theta_chain_ibp <- s_chain_ibp - alpha_chain_ibp
    
    params_ibp[[lab_a]] <- a_chain_ibp
    params_ibp[[lab_b]] <- b_chain_ibp
    params_ibp[[lab_alpha]] <- alpha_chain_ibp
    params_ibp[[lab_theta]] <- theta_chain_ibp
    
    # ############ 5.1) Run SB-SP #################
    # 
    # # Set tau for MALA
    # tau_sp <- 0.0001
    # 
    # output_sp <- gibbs_sampler_sb_sp(Z = train_mat,
    #                                  c_0_sp, beta_0_sp, alpha_0_sp,
    #                                  p_sp, r_sp, t_sp, a_alpha_sp, b_alpha_sp,
    #                                  tau_sp, fixed = c(F,F,F),
    #                                  S_ibp, n_burnin_ibp, thin_ibp, seed)
    # 
    # n_saved_iter_sp <- length(output_sp$c_vec)
    # c_chain_sp <- output_sp$c_vec
    # beta_chain_sp <- output_sp$beta_vec
    # alpha_chain_sp <- output_sp$alpha_vec
    # 
    # params_sp[[paste0("c.",N)]] <- c_chain_sp
    # params_sp[[paste0("beta.",N)]] <- beta_chain_sp
    # params_sp[[paste0("alpha.",N)]] <- alpha_chain_sp
    
    ####### 6) Estimate limit distributions (Poiss/NB) ################
    ntilde_chain_poiss <- generate_Ntilde_chain_poiss(lambda_poiss, alpha_chain_poiss,
                                                      theta_chain_poiss, n = N, Kn)
    
    ntilde_chain_negbin <- generate_Ntilde_chain_negbin(nstar_nb, p_nb,
                                                        alpha_chain_negbin,
                                                        theta_chain_negbin, n = N, Kn)
    
    # ntilde_chain_negbin_prior <- generate_Ntilde_chain_negbin(nstar_chain_negbin_prior, 
    #                                                           p_chain_negbin_prior,
    #                                                           alpha_chain_negbin,
    #                                                           theta_chain_negbin, n = N, Kn)
    
    gg_ntilde_poiss[[lab_comb]] <- ntilde_chain_poiss
    gg_ntilde_negbin[[lab_comb]] <- ntilde_chain_negbin
    #gg_ntilde_negbin_prior[[paste0("N.",N)]] <- ntilde_chain_negbin_prior
    
    ######## 7) Extrapolation in the test set (Poiss/NB/Gamma/SP) ##############
    # Poisson
    kmn_chain_poiss <- generate_Kmn_chain_poiss(lambda_poiss, alpha_chain_poiss,
                                                theta_chain_poiss, M = M, n = N)
    
    est_ci_pred_poiss <- matrix(NA, nrow = M, ncol = 3)
    # first column = lower bound
    # second columns = medians
    # third columns = upper bound
    for (m in 1:M){
      est_ci_pred_poiss[m,] <- quantile(kmn_chain_poiss[m,], probs = c(0.025,0.5,0.975))
    }
    est_ci_pred_poiss <- list("medians" = est_ci_pred_poiss[,2],
                              "lbs" = est_ci_pred_poiss[,1],
                              "ubs" = est_ci_pred_poiss[,3])
    
    list_kmn_pred_test_poiss[[lab_comb]] <- est_ci_pred_poiss
    
    # Negative Binomial (fixed params)
    kmn_chain_negbin <- generate_Kmn_chain_negbin(nstar_nb, p_nb,
                                                  alpha_chain_negbin, theta_chain_negbin,
                                                  M = M, n = N, Kn)
    
    est_ci_pred_negbin <- matrix(NA, nrow = M, ncol = 3)
    # first column = lower bound
    # second columns = medians
    # third columns = upper bound
    for (m in 1:M){
      est_ci_pred_negbin[m,] <- quantile(kmn_chain_negbin[m,], probs = c(0.025,0.5,0.975))
    }
    est_ci_pred_negbin <- list("medians" = est_ci_pred_negbin[,2],
                               "lbs" = est_ci_pred_negbin[,1],
                               "ubs" = est_ci_pred_negbin[,3])
    
    list_kmn_pred_test_negbin[[lab_comb]] <- est_ci_pred_negbin
    
    # # Negative Binomial (prior on params)
    # kmn_chain_negbin_prior <- generate_Kmn_chain_negbin(nstar_chain_negbin_prior,
    #                                                     p_chain_negbin_prior,
    #                                                     alpha_chain_negbin, 
    #                                                     theta_chain_negbin,
    #                                                     M = M, n = N, Kn)
    # 
    # est_ci_pred_negbin_prior <- matrix(NA, nrow = M, ncol = 3)
    # # first column = lower bound
    # # second columns = medians
    # # third columns = upper bound
    # for (m in 1:M){
    #   est_ci_pred_negbin_prior[m,] <- quantile(kmn_chain_negbin_prior[m,], probs = c(0.025,0.5,0.975))
    # }
    # est_ci_pred_negbin_prior <- list("medians" = est_ci_pred_negbin_prior[,2],
    #                            "lbs" = est_ci_pred_negbin_prior[,1],
    #                            "ubs" = est_ci_pred_negbin_prior[,3])
    # 
    # list_kmn_pred_test_negbin_prior[[paste0("N.",N)]] <- est_ci_pred_negbin_prior
    
    
    # IBP + Gamma
    kmn_chain_ibp <- generate_Kmn_chain_gamma_ibp(a_chain_ibp, b_chain_ibp, alpha_chain_ibp,
                                                  theta_chain_ibp, M = M, n = N, Kn)
    
    est_ci_pred_ibp <- matrix(NA, nrow = M, ncol = 3)
    # first column = lower bound
    # second columns = medians
    # third columns = upper bound
    for (m in 1:M){
      est_ci_pred_ibp[m,] <- quantile(kmn_chain_ibp[m,], probs = c(0.025,0.5,0.975))
    }
    est_ci_pred_ibp <- list("medians" = est_ci_pred_ibp[,2],
                            "lbs" = est_ci_pred_ibp[,1],
                            "ubs" = est_ci_pred_ibp[,3])
    
    list_kmn_pred_test_ibp[[lab_comb]] <- est_ci_pred_ibp
    
    # # SB-SP
    # kmn_chain_sp <- generate_Kmn_chain_gamma_ibp(a_chain = c_chain_sp + 1,
    #                                              b_chain = beta_chain_sp*(1-alpha_chain_sp)/alpha_chain_sp, 
    #                                              alpha_chain = alpha_chain_sp,
    #                                              theta_chain = 1 - alpha_chain_sp,
    #                                              M = M, n = N, Kn)
    # 
    # est_ci_pred_sp <- matrix(NA, nrow = M, ncol = 3)
    # # first column = lower bound
    # # second columns = medians
    # # third columns = upper bound
    # for (m in 1:M){
    #   est_ci_pred_sp[m,] <- quantile(kmn_chain_sp[m,], probs = c(0.025,0.5,0.975))
    # }
    # est_ci_pred_sp <- list("medians" = est_ci_pred_sp[,2],
    #                         "lbs" = est_ci_pred_sp[,1],
    #                         "ubs" = est_ci_pred_sp[,3])
    # 
    # list_kmn_pred_test_sp[[paste0("N.",N)]] <- est_ci_pred_sp
    
    
  }
  
}



############# 8) Save results:  MCMC convergence ####################

save(params_poiss, file = "m4_params_poiss.Rda")
save(params_negbin, file = "m4_params_negbin.Rda")
#save(params_negbin_prior, file = "m4_params_negbin_prior.Rda")
save(params_ibp, file = "m4_params_ibp.Rda")
#save(params_sp, file = "m4_params_sp.Rda")

############ 9) Save results: samples from limiting distributions (Poiss/NB) ##############
# Poisson
save(gg_ntilde_poiss, file = "m4_ntilde_poiss.Rda")
# Negative Binomial (fixed)
save(gg_ntilde_negbin, file = "m4_ntilde_negbin.Rda")
# Negative Binomial (prior)
#save(gg_ntilde_negbin_prior, file = "m4_ntilde_negbin_prior.Rda")

######### 10) Save results: CI for extrapolation (Poiss/NB/Gamma/SP) ################
# Poisson
saveRDS(list_kmn_pred_test_poiss, "m4_ci_poiss.rds")
# Negative Binomial (fixed)
saveRDS(list_kmn_pred_test_negbin, "m4_ci_negbin.rds")
# Negative Binomial (prior)
#saveRDS(list_kmn_pred_test_negbin_prior, "m4_ci_negbin_prior.rds")
# Gamma IBP
saveRDS(list_kmn_pred_test_ibp, "m4_ci_ibp.rds")
# SB-SP
#saveRDS(list_kmn_pred_test_sp, "m4_ci_sp.rds")


######## 11) Save the data ##############################
saveRDS(data_mat, "m4_data_mat.rds")





# ##### Accuracy on multiple datasets #####
# 
# # number of datasets to average over
# D <- 50
# 
# # mcmc parameters
# S_poiss <- S_negbin <- S_ibp <- 5*10^4
# n_burnin_poiss <- n_burnin_negbin <- n_burnin_ibp<- 10^4
# thin_poiss <- thin_negbin <- thin_ibp <- 2
# number_saved_iterations_poiss <- (S_poiss - n_burnin_poiss)/thin_poiss 
# number_saved_iterations_negbin <- (S_negbin - n_burnin_negbin)/thin_negbin 
# number_saved_iterations_ibp <- (S_ibp - n_burnin_ibp)/thin_ibp 
# 
# # store objects
# avg_ntilde_poiss <- data.frame(matrix(nrow = D, ncol = length(Ns)))
# colnames(avg_ntilde_poiss) <- paste("N", Ns, sep = ".")
# avg_ntilde_negbin <- data.frame(matrix(nrow = D, ncol = length(Ns)))
# colnames(avg_ntilde_negbin) <- paste("N", Ns, sep = ".")
# avg_ntilde_negbin_prior <- data.frame(matrix(nrow = D, ncol = length(Ns)))
# colnames(avg_ntilde_negbin_prior) <- paste("N", Ns, sep = ".")
# 
# # number of new features observed in the test
# obs_new <- data.frame(matrix(nrow = D, ncol = length(Ns)))
# colnames(obs_new) <- paste("N", Ns, sep = ".")
# 
# # number of old features observed in the training
# obs_train <- data.frame(matrix(nrow = D, ncol = length(Ns)))
# colnames(obs_train) <- paste("N", Ns, sep = ".")
# 
# est_new_poiss <- data.frame(matrix(nrow = D, ncol = length(Ns)))
# colnames(est_new_poiss) <- paste("N", Ns, sep = ".")
# est_new_negbin <- data.frame(matrix(nrow = D, ncol = length(Ns)))
# colnames(est_new_negbin) <- paste("N", Ns, sep = ".")
# est_new_negbin_prior <- data.frame(matrix(nrow = D, ncol = length(Ns)))
# colnames(est_new_negbin_prior) <- paste("N", Ns, sep = ".")
# est_new_ibp <- data.frame(matrix(nrow = D, ncol = length(Ns)))
# colnames(est_new_ibp) <- paste("N", Ns, sep = ".")
# est_new_sp <- data.frame(matrix(nrow = D, ncol = length(Ns)))
# colnames(est_new_sp) <- paste("N", Ns, sep = ".")
# 
# seed = 123456
# set.seed(seed)
# 
# for (d in 1:D){
#   data_mat <- matrix(rbinom(L*H, size = 1, prob = rep(pis, L)), 
#                      nrow = L, ncol = H, byrow = T )
#   data_list <- create_features_list(data_mat)
#   
#   for (j in 1:length(Ns)){
#     N <- Ns[j]
#     M <- L - N
#     
#     train_mat <- data_mat[1:N,]
#     test_mat <- data_mat[(N+1):L, ]
#     # convert the binary matrix into list of features
#     train_list <- create_features_list(train_mat)
#     test_list <- create_features_list(test_mat)
#     
#     ################# 3) Run BB + Poisson ###########
#     
#     # Set tau for MALA
#     tau_poiss <- 0.0005
#     
#     output_poiss <- gibbs_sampler_poiss_fixed_lambda(Z = train_mat,
#                                                      alpha_bar_0_bb, s_0_bb,
#                                                      lambda_poiss, a_alpha_bb, b_alpha_bb, a_s_bb, b_s_bb,
#                                                      tau_poiss,
#                                                      S_poiss, n_burnin_poiss, thin_poiss, seed)
#     
#     n_saved_iter_poiss <- length(output_poiss$s_vec)
#     s_chain_poiss <- output_poiss$s_vec
#     alpha_bar_chain_poiss <- output_poiss$alpha_bar_vec
#     alpha_chain_poiss <- - alpha_bar_chain_poiss
#     theta_chain_poiss <- s_chain_poiss + alpha_bar_chain_poiss
#     
#     
#     ####### 4) Run BB + Negative-Binomial (fixed parameters) ############
#     
#     # Set tau for MALA
#     tau_nb <- 0.001
#     
#     output_negbin <- gibbs_sampler_negbin_geometric_fixed_pars(Z = train_mat,
#                                                                s_0_bb, alpha_bar_0_bb,
#                                                                nstar_nb, p_nb,
#                                                                a_alpha_bb, b_alpha_bb, a_s_bb, b_s_bb,
#                                                                tau_nb, 
#                                                                S_negbin, n_burnin_negbin, thin_negbin, seed)
#     
#     n_saved_iter_negbin <- length(output_negbin$s_vec)
#     s_chain_negbin <- output_negbin$s_vec
#     alpha_bar_chain_negbin <- output_negbin$alpha_bar_vec
#     alpha_chain_negbin <- - alpha_bar_chain_negbin
#     theta_chain_negbin <- s_chain_negbin + alpha_bar_chain_negbin
#     
#     
#     ####### 4.b) Run BB + Negative-Binomial (prior on parameters) #######
#     
#     # Set tau for MALA
#     tau_nb_prior <- 0.001
#     
#     output_negbin_prior <- gibbs_sampler_negbin_geometric_prior_pars(Z = train_mat,
#                                                                      nstar_0_nb, p_0_nb, s_0_bb, alpha_bar_0_bb,
#                                                                      q_star_nb, alpha_p_nb, beta_p_nb,
#                                                                      a_alpha_bb, b_alpha_bb, a_s_bb, b_s_bb,
#                                                                      tau_nb_prior, fixed = c(F,F,F,F),
#                                                                      S_negbin, n_burnin_negbin, thin_negbin, seed)
#     
#     n_saved_iter_negbin_prior <- length(output_negbin_prior$nstar_vec)
#     nstar_chain_negbin_prior <- output_negbin_prior$nstar_vec
#     p_chain_negbin_prior <- output_negbin_prior$p_vec
#     s_chain_negbin_prior <- output_negbin_prior$s_vec
#     alpha_bar_chain_negbin_prior <- output_negbin_prior$alpha_bar_vec
#     alpha_chain_negbin_prior <- - alpha_bar_chain_negbin_prior
#     theta_chain_negbin_prior <- s_chain_negbin_prior + alpha_bar_chain_negbin_prior
#     
#     
#     ############ 5) Run IBP + Gamma #################
#     
#     # Set tau for MALA
#     tau_ibp <- 0.002
#     
#     output_ibp <- gibbs_sampler_gamma_ibp(Z = train_mat,
#                                           a_0_ibp, b_0_ibp, s_0_ibp, alpha_0_ibp,
#                                           p_ibp, r_ibp, t_ibp, a_alpha_ibp, b_alpha_ibp, a_s_ibp, b_s_ibp,
#                                           tau_ibp, fixed = c(F,F,F,F),
#                                           S_ibp, n_burnin_ibp, thin_ibp, seed)
#     
#     n_saved_iter_ibp <- length(output_ibp$a_vec)
#     a_chain_ibp <- output_ibp$a_vec
#     b_chain_ibp <- output_ibp$b_vec
#     s_chain_ibp <- output_ibp$s_vec
#     alpha_chain_ibp <- output_ibp$alpha_vec
#     theta_chain_ibp <- s_chain_ibp - alpha_chain_ibp
#     
#     ############ 5.1) Run SB-SP #################
#     
#     # Set tau for MALA
#     tau_sp <- 0.01
#     
#     output_sp <- gibbs_sampler_sb_sp(Z = train_mat,
#                                      c_0_sp, beta_0_sp, alpha_0_sp,
#                                      p_sp, r_sp, t_sp, a_alpha_sp, b_alpha_sp,
#                                      tau_sp, fixed = c(F,F,F),
#                                      S_ibp, n_burnin_ibp, thin_ibp, seed)
#     
#     n_saved_iter_sp <- length(output_sp$c_vec)
#     c_chain_sp <- output_sp$c_vec
#     beta_chain_sp <- output_sp$beta_vec
#     alpha_chain_sp <- output_sp$alpha_vec
#     
#     
#     ##### 6) Estimate mean of the limit distributions (Poiss/NB) ################
#     Kn = ncol(train_mat[,colSums(train_mat) > 0])
#     
#     ntilde_chain_poiss <- generate_Ntilde_chain_poiss(lambda_poiss, alpha_chain_poiss,
#                                                       theta_chain_poiss, n = N, Kn)
#     
#     ntilde_chain_negbin <- generate_Ntilde_chain_negbin(nstar_nb, p_nb,
#                                                         alpha_chain_negbin,
#                                                         theta_chain_negbin, n = N, Kn)
#     
#     ntilde_chain_negbin_prior <- generate_Ntilde_chain_negbin(nstar_chain_negbin_prior, p_chain_negbin_prior,
#                                                               alpha_chain_negbin,
#                                                               theta_chain_negbin, n = N, Kn)
#     
#     avg_ntilde_poiss[d, j] <- mean(ntilde_chain_poiss)
#     avg_ntilde_negbin[d, j] <- mean(ntilde_chain_negbin)
#     avg_ntilde_negbin_prior[d, j] <- mean(ntilde_chain_negbin_prior)
#     
#     #### 7) Collect quantities for accuracy (Poiss/NB/Gamma/SP) #####
#     feat_train <- unique(unlist(train_list))
#     feat_test <- unique(unlist(test_list))
#     # number of features observed in training
#     obs_train[d,j] <- length(feat_train)
#     # number of new features observed in test
#     obs_new_features <- setdiff(feat_test, feat_train)
#     obs_new[d,j] <- length(obs_new_features)
#     
#     # Poisson
#     kmn_chain_poiss <- generate_Kmn_chain_poiss(lambda_poiss, alpha_chain_poiss,
#                                                 theta_chain_poiss, M = M, n = N)
#     
#     est_new_poiss[d , j] <- mean(kmn_chain_poiss)
#     
#     # Negative Binomial (fixed)
#     kmn_chain_negbin <- generate_Kmn_chain_negbin(nstar_nb, p_nb,
#                                                   alpha_chain_negbin, theta_chain_negbin,
#                                                   M = M, n = N, Kn)
#     
#     est_new_negbin[d , j] <- mean(kmn_chain_negbin)
#     
#     # Negative Binomial (prior)
#     kmn_chain_negbin_prior <- generate_Kmn_chain_negbin(nstar_chain_negbin_prior, p_chain_negbin_prior,
#                                                         alpha_chain_negbin, theta_chain_negbin,
#                                                         M = M, n = N, Kn)
#     
#     est_new_negbin_prior[d , j] <- mean(kmn_chain_negbin_prior)
#     
#     # IBP + Gamma
#     kmn_chain_ibp <- generate_Kmn_chain_gamma_ibp(a_chain_ibp, b_chain_ibp, alpha_chain_ibp,
#                                                   theta_chain_ibp, M = M, n = N, Kn)
#     
#     est_new_ibp[d , j] <- mean(kmn_chain_ibp)
#     
#     # SB-SP
#     kmn_chain_sp <- generate_Kmn_chain_gamma_ibp(a_chain = c_chain_sp + 1,
#                                                  b_chain = beta_chain_sp*(1-alpha_chain_sp)/alpha_chain_sp, 
#                                                  alpha_chain = alpha_chain_sp,
#                                                  theta_chain = 1 - alpha_chain_sp,
#                                                  M = M, n = N, Kn)
#     
#     est_new_sp[d , j] <- mean(kmn_chain_sp)
#     
#   }
#   
# }
# 
# 
# ###### 8) Save results: limit distribution estimates #####
# save(avg_ntilde_poiss, file = "m4_avg_ntilde_poiss.Rda")
# save(avg_ntilde_negbin, file = "m4_avg_ntilde_negbin.Rda")
# save(avg_ntilde_negbin_prior, file = "m4_avg_ntilde_negbin_prior.Rda")
# 
# 
# ###### 9) Save results: quantities for accuracy #####
# save(obs_train, file = "m4_obs_train.Rda")
# save(obs_new, file = "m4_obs_new.Rda")
# save(est_new_poiss, file = "m4_est_new_poiss.Rda")
# save(est_new_negbin, file = "m4_est_new_negbin.Rda")
# save(est_new_negbin_prior, file = "m4_est_new_negbin_prior.Rda")
# save(est_new_ibp, file = "m4_est_new_ibp.Rda")
# save(est_new_sp, file = "m4_est_new_sp.Rda")
# 
# 
