
#' Function the gradient of the log full-conditional for s_hat and alpha_bar_hat
#' 
#' @param nstar_hat
#' @param s_hat 
#' @param alpha_bar_hat 
#' @param lambda 
#' @param n 
#' @param K 
#' @param counts 
#' @param a_s 
#' @param b_s 
#' @param a_alpha 
#' @param b_alpha 
#'
#' @export
#'
compute_grad_log_full_negbin <- function(nstar_hat, s_hat, alpha_bar_hat, p,
                                         n, K, counts, a_star, b_star,
                                         a_s, b_s, a_alpha, b_alpha){
  
  enstar <- exp(nstar_hat)
  es <- exp(s_hat)
  ea_bar <- exp(alpha_bar_hat)
  
  p1 <- exp(lgamma(n+es) + lgamma(es + ea_bar) -
              lgamma(es) - lgamma(es + ea_bar +n) )
  
  # derivative wrt nstar_hat
  dnstar_hat <- a_star + enstar * (digamma(K + enstar) - digamma(enstar) - 
                                     log(1-(1-p)*p1) + log(p) - b_star)
    
  # derivative wrt s_hat
  
  p2 <- exp( log(enstar + K) + log(1-p) - log(1/p1 -1 + p) )
  
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
  
  return (c(dnstar_hat, ds_hat, dalpha_bar_hat))
}



############################################################################


#' Function computing the log ratio of the full-conditional of s_hat, alpha_bar_hat
#' in proposed and current values
#'
#' @param s_hat_prop 
#' @param alpha_bar_hat_prop 
#' @param s_hat_curr 
#' @param alpha_bar_hat_curr 
#' @param lambda 
#' @param n 
#' @param K 
#' @param counts 
#' @param a_s 
#' @param b_s 
#' @param a_alpha 
#' @param b_alpha 
#'
#' @export
#'
compute_log_ratio_full_negbin <- function(nstar_hat_prop, s_hat_prop, alpha_bar_hat_prop,
                                          nstar_hat_curr, s_hat_curr, alpha_bar_hat_curr, 
                                          p, n, K, counts, a_star, b_star, 
                                          a_s, b_s, a_alpha, b_alpha){
  
  enstar_prop <- exp(nstar_hat_prop)
  es_prop <- exp(s_hat_prop)
  ea_bar_prop <- exp(alpha_bar_hat_prop)
  enstar_curr <- exp(nstar_hat_curr)
  es_curr <- exp(s_hat_curr)
  ea_bar_curr <- exp(alpha_bar_hat_curr)
  
  p1_prop <- exp(lgamma(n+es_prop) + lgamma(es_prop + ea_bar_prop) -
                   lgamma(es_prop) - lgamma(es_prop + ea_bar_prop +n) )
  p1_curr <- exp(lgamma(n+es_curr) + lgamma(es_curr + ea_bar_curr) -
                   lgamma(es_curr) - lgamma(es_curr + ea_bar_curr +n) )
  
  v <- lgamma(ea_bar_prop + counts) - lgamma(1 + ea_bar_prop) + lgamma(es_prop+n - counts) -
    lgamma(es_prop) - lgamma(ea_bar_curr + counts) + lgamma(1 + ea_bar_curr) - 
    lgamma(es_curr+n - counts) + lgamma(es_curr)
  
  res <- lgamma(K + enstar_prop) - lgamma(K+enstar_curr) - 
    lgamma(enstar_prop) + lgamma(enstar_prop) + 
    K*(lgamma(n+es_curr+ea_bar_curr) + lgamma(es_prop + ea_bar_prop) -
         lgamma(es_curr + ea_bar_curr) - lgamma(es_prop+ea_bar_prop+n)) -
    (enstar_prop + K)*log(1-(1-p)*p1_prop) + (enstar_curr + K)*log(1-(1-p)*p1_curr) +
    sum(v) + 
    (nstar_hat_prop - nstar_hat_curr)*a_star + (enstar_prop - enstar_curr)*(log(p) - b_star) + 
    (K+a_alpha)*(alpha_bar_hat_prop - alpha_bar_hat_curr) - b_alpha*(ea_bar_prop - ea_bar_curr) + 
    a_s*(s_hat_prop - s_hat_curr) - b_s*(es_prop - es_curr) 
    
    
  return(res)
  
}


#############################################################################

#' Function computing the log of the ratio of the terms related to q in MALA
#'
#' @param s_hat_prop 
#' @param alpha_bar_hat_prop 
#' @param s_hat_curr 
#' @param alpha_bar_hat_curr 
#' @param tau 
#' @param grad_log_full_curr 
#' @param grad_log_full_prop 
#'
#' @export
#'
compute_log_ratio_q_negbin <- function(nstar_hat_prop, s_hat_prop, alpha_bar_hat_prop,
                                       nstar_hat_curr, s_hat_curr, alpha_bar_hat_curr, 
                                       tau, grad_log_full_curr, grad_log_full_prop){
  
  norm_prop <- norm(c(nstar_hat_curr, s_hat_curr, alpha_bar_hat_curr) - 
                      c(nstar_hat_prop, s_hat_prop, alpha_bar_hat_prop) -
                      tau*grad_log_full_prop, type="2")^2
  
  norm_curr <- norm( c(nstar_hat_prop, s_hat_prop, alpha_bar_hat_prop) - 
                       c(nstar_hat_curr, s_hat_curr, alpha_bar_hat_curr) -
                       tau*grad_log_full_curr, type="2")^2
  
  res <- 1/(4*tau)*(norm_curr - norm_prop)
  
  return (res)
}


############################################################################


#' Metropolis-within-Gibbs sampler for BB + NegBin (Gamma prior on nstar)
#'
#' @param Z [integer] binary matrix of presence/absence (n x K- dimensional)
#' @param nstar_0 [numeric] initial value of nstar
#' @param p_0 [numeric] initial value of p
#' @param s_0 [numeric] initial value of s
#' @param alpha_bar_0 [numeric] initial value of alpha_bar 
#' @param a_star [numeric] 
#' @param b_star [numeric]
#' @param alpha_p [numeric] 
#' @param beta_p [numeric]
#' @param a_alpha [numeric]
#' @param b_alpha [numeric]
#' @param a_s [numeric]
#' @param b_s [numeric]
#' @param tau [numeric] MALA step-size
#' @param fixed [logical] vector indicating if the variable is fixed or has prior:
#' if fixed, it stays equal to initial value
#' @param S [integer] number of iterations for the MCMC algorithm
#' @param n_burnin [integer] number of iterations for the burn-in
#' @param thin [integer] thinning
#' @param seed [integer] seed
#'
#' @return
#' @export
#'
#' @examples
gibbs_sampler_negbin_gamma <- function(Z, 
                                nstar_0, p_0, s_0, alpha_bar_0,
                                a_star, b_star, alpha_p, beta_p, a_alpha, b_alpha, a_s, b_s,
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
    
    
    ################################################################
    ############# Draw p | Z, alpha_bar, s, nstar ###################
    ###############################################################
    
    if (fixed[2] == FALSE){
      # Update Ntilde | p, Z, alpha_bar, s, nstar
      t <- lgamma(s + n) + lgamma(s + alpha_bar) - lgamma(s) - lgamma(s+alpha_bar+n)
      
      Ntilde <- K + rnbinom(1, size = nstar + K, prob = 1 - (1-p)* exp(t))
      
      # Update p | Ntilde, Z, alpha_bar, s, nstar
      
      p <- rbeta(1, shape1 = nstar + alpha_p, shape2 = Ntilde + beta_p)
      
    }
    
    ################################################################
    ############# Draw (nstar, s, alpha_bar) | Z, p ##################
    ###############################################################
    
    # In order to update nstar, s, alpha_bar, we update nstar_hat, s_hat, alpha_bar_hat, 
    # defined as the logarithm of nstar, s, alpha_bar
    
    ### Current values for nstar_hat, s_hat, alpha_bar_hat
    nstar_hat_curr <- log(nstar)
    s_hat_curr <- log(s)
    alpha_bar_hat_curr <- log(alpha_bar)
    
    ### Propose values for nstar_hat, s_hat, alpha_bar_hat
    # Compute the gradient of log-full conditional for MALA 
    # (the order of returned variables is: nstar_hat, s_hat, alpha_bar_hat)
    grad_log_full_curr <- compute_grad_log_full_negbin(nstar_hat_curr, s_hat_curr, alpha_bar_hat_curr, 
                                                       p, n, K, counts, a_star, b_star,
                                                       a_s, b_s, a_alpha, b_alpha)
    if (fixed[1] == TRUE){
      nstar_hat_prop <- nstar_hat_curr
    } else {
      nstar_hat_prop <- nstar_hat_curr + tau*grad_log_full_curr[1] + sqrt(2*tau)*rnorm(1)
    }
    
    if (fixed[3] == TRUE){
      s_hat_prop <- s_hat_curr
    } else {
      s_hat_prop <- s_hat_curr + tau*grad_log_full_curr[2] + sqrt(2*tau)*rnorm(1)
    }
    
    if (fixed[4] == TRUE){
      alpha_bar_hat_prop <- alpha_bar_hat_curr
    } else {
      alpha_bar_hat_prop <- alpha_bar_hat_curr + tau*grad_log_full_curr[3] + sqrt(2*tau)*rnorm(1)
    }
    
    ### Acceptance probability 
    # Compute the log ratio of the full-cond in prop point and curr point
    log_ratio_full <- compute_log_ratio_full_negbin(nstar_hat_prop, s_hat_prop, alpha_bar_hat_prop,
                                                   nstar_hat_curr, s_hat_curr, alpha_bar_hat_curr, 
                                                   p, n, K, counts, a_star, b_star, 
                                                   a_s, b_s, a_alpha, b_alpha)
    
    # Compute the log ratio of the terms related to the proposal q
    grad_log_full_prop <- compute_grad_log_full_negbin(nstar_hat_prop, s_hat_prop, alpha_bar_hat_prop, 
                                                       p, n, K, counts, a_star, b_star,
                                                       a_s, b_s, a_alpha, b_alpha)
    
    log_ratio_q <- compute_log_ratio_q_negbin(nstar_hat_prop, s_hat_prop, alpha_bar_hat_prop,
                                              nstar_hat_curr, s_hat_curr, alpha_bar_hat_curr, tau,
                                              grad_log_full_curr, grad_log_full_prop)
    
    # Compute acceptance probability
    log_acc_prob <- log_ratio_full + log_ratio_q
    acc_prob <- min(1, exp(log_acc_prob))
    
    # Decide if accept or not the new parameter vector
    if (runif(1) < acc_prob){ # accept
      nstar <- exp(nstar_hat_prop)
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


############################################################################
############################################################################
############################################################################

#' Function the gradient of the log full-conditional for s_hat and alpha_bar_hat
#' 
#' @param nstar_hat
#' @param s_hat 
#' @param alpha_bar_hat 
#' @param lambda 
#' @param n 
#' @param K 
#' @param counts 
#' @param a_s 
#' @param b_s 
#' @param a_alpha 
#' @param b_alpha 
#'
#' @export
#'
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

############################################################################


#' Function computing the log ratio of the full-conditional of s_hat, alpha_bar_hat
#' in proposed and current values
#'
#' @param s_hat_prop 
#' @param alpha_bar_hat_prop 
#' @param s_hat_curr 
#' @param alpha_bar_hat_curr 
#' @param lambda 
#' @param n 
#' @param K 
#' @param counts 
#' @param a_s 
#' @param b_s 
#' @param a_alpha 
#' @param b_alpha 
#'
#' @export
#'
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


#############################################################################

#' Function computing the log of the ratio of the terms related to q in MALA
#'
#' @param s_hat_prop 
#' @param alpha_bar_hat_prop 
#' @param s_hat_curr 
#' @param alpha_bar_hat_curr 
#' @param tau 
#' @param grad_log_full_curr 
#' @param grad_log_full_prop 
#'
#' @export
#'
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

######################################################################


#' Metropolis-within-Gibbs sampler for BB + NegBin (Geometric prior on nstar)
#'
#' @param Z [integer] binary matrix of presence/absence (n x K- dimensional)
#' @param nstar_0 [numeric] initial value of nstar
#' @param p_0 [numeric] initial value of p
#' @param s_0 [numeric] initial value of s
#' @param alpha_bar_0 [numeric] initial value of alpha_bar 
#' @param q_star [numeric] hyperparameter of Geometric prior on nstar
#' @param alpha_p [numeric] 
#' @param beta_p [numeric]
#' @param a_alpha [numeric]
#' @param b_alpha [numeric]
#' @param a_s [numeric]
#' @param b_s [numeric]
#' @param tau [numeric] MALA step-size
#' @param fixed [logical] vector indicating if the variable is fixed or has prior:
#' if fixed, it stays equal to initial value
#' @param S [integer] number of iterations for the MCMC algorithm
#' @param n_burnin [integer] number of iterations for the burn-in
#' @param thin [integer] thinning
#' @param seed [integer] seed
#'
#' @return
#' @export
#'
#' @examples
gibbs_sampler_negbin_geometric <- function(Z, 
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


#########################################################################
#########################################################################



#' Generate the chain for Kmn, from m=1 to m=M, given the output chains of the mcmc
#'
#' @param nstar_chain
#' @param p_chain
#' @param alpha_chain
#' @param theta_chain
#' @param M
#' @param n
#' @param Kn
#'
#' @export
#'
generate_Kmn_chain_negbin <- function(nstar_chain, p_chain, alpha_chain, theta_chain, M, n, Kn = 0){
  
  if (n == 0 & Kn != 0){
    stop("if n=0, the number of observed features Kn must be 0!")
  }
  
  S <- length(nstar_chain)
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




##########################################

#' Title
#'
#' @param lambda_chain_poiss 
#' @param alpha_chain_poiss 
#' @param theta_chain_poiss 
#' @param n 
#' @param Kn
#'
#' @export
#'
generate_Ntilde_chain_negbin <- function(nstar_chain_negbin, p_chain_negbin,
                                         alpha_chain_negbin, theta_chain_negbin,
                                         n, Kn){
  
  S <- length(nstar_chain_negbin)
  
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
