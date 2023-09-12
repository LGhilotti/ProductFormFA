
############################################################################
############################################################################
############################################################################

#' Function for the log posterior density of alpha_hat and s_hat
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

############################################################################

#' Function the gradient of the log full-conditional for s_hat and alpha_bar_hat
#' 
#' @param n 
#' @param K 
#' @param counts 
#' @param a_alpha 
#' @param b_alpha 
#'
#' @export
#'
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


#############################################################################

#' Function computing the log ratio of the full-conditional of s_hat, alpha_bar_hat
#' in proposed and current values
#'
#' @param s_hat_prop 
#' @param alpha_hat_prop 
#' @param s_hat_curr 
#' @param alpha_hat_curr 
#' @param Gam 
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


#############################################################################

#' Function computing the log ratio of the full-conditional of s_hat, alpha_bar_hat
#' in proposed and current values
#'
#' @param Gam 
#' @param n 
#' @param K 
#' @param counts 
#' @param a_alpha 
#' @param b_alpha 
#'
#' @export
#'
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
#############################################################################

#' Function computing the log of the ratio of the terms related to q in MALA
#'
#' @param s_hat_prop 
#' @param alpha_hat_prop 
#' @param s_hat_curr 
#' @param alpha_hat_curr 
#' @param tau 
#' @param grad_log_full_curr 
#' @param grad_log_full_prop 
#'
#' @export
#'
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



######################################################################

#' Function computing the log of the ratio of the terms related to q in MALA
#'
#' @param s_hat_prop 
#' @param alpha_hat_prop 
#' @param s_hat_curr 
#' @param alpha_hat_curr 
#' @param tau 
#' @param grad_log_full_curr 
#' @param grad_log_full_prop 
#'
#' @export
#'
compute_log_ratio_q_sb_sp <- function(alpha_hat_prop,
                                      alpha_hat_curr, tau,
                                      der_log_full_curr, der_log_full_prop){
  
  norm_prop <- (alpha_hat_curr - alpha_hat_prop - tau*der_log_full_prop)^2
  
  norm_curr <- (alpha_hat_prop - alpha_hat_curr - tau*der_log_full_curr)^2
  
  res <- 1/(4*tau)*(norm_curr - norm_prop)
  
  return (res)
  
}

######################################################################


#' Metropolis-within-Gibbs sampler for BB + NegBin (Geometric prior on nstar)
#'
#' @param Z [integer] binary matrix of presence/absence (n x K- dimensional)
#' @param a_0 [numeric] initial value of a
#' @param b_0 [numeric] initial value of b
#' @param s_0 [numeric] initial value of s
#' @param alpha_0 [numeric] initial value of alpha 
#' @param p [numeric] hyperparameter of Geometric prior on a
#' @param r [numeric] 
#' @param t [numeric]
#' @param a_alpha [numeric]
#' @param b_alpha [numeric]
#' @param a_s [numeric]
#' @param b_s [numeric]
#' @param tau [numeric] MALA step-size
#' @param fixed [logical] vector indicating if the variable is fixed or has prior:
#' if fixed, it stays equal to initial value - order: a, b, alpha, s
#' @param S [integer] number of iterations for the MCMC algorithm
#' @param n_burnin [integer] number of iterations for the burn-in
#' @param thin [integer] thinning
#' @param seed [integer] seed
#'
#' @return
#' @import numDeriv
#' @import stats
#' @import MASS
#' @import matlib
#' @export
#'
gibbs_sampler_gamma_ibp_mala <- function(Z, 
                                    a_0, b_0, s_0, alpha_0,
                                    p, r, t, a_alpha, b_alpha, a_s, b_s,
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
    
    ### Propose values for s_hat, alpha_hat
    # Compute the gradient of log-full conditional for MALA 
    # (the order of returned variables is: s_hat, alpha_hat)
    grad_log_full_curr <- compute_grad_log_full_gamma_ibp(s_hat_curr, alpha_hat_curr, 
                                                            Gam, n, K, counts, 
                                                            a_s, b_s, a_alpha, b_alpha)
    
    if (fixed[4] == TRUE){
      s_hat_prop <- s_hat_curr
    } else {
      s_hat_prop <- s_hat_curr + tau*grad_log_full_curr[1] + sqrt(2*tau)*rnorm(1)
    }
    
    if (fixed[3] == TRUE){
      alpha_hat_prop <- alpha_hat_curr
    } else {
      alpha_hat_prop <- alpha_hat_curr + tau*grad_log_full_curr[2] + sqrt(2*tau)*rnorm(1)
    }
    
    ### Acceptance probability 
    # Compute the log ratio of the full-cond in prop point and curr point
    log_ratio_full <- compute_log_ratio_full_gamma_ibp(s_hat_prop, alpha_hat_prop,
                                                         s_hat_curr, alpha_hat_curr, 
                                                         Gam, n, K, counts,   
                                                         a_s, b_s, a_alpha, b_alpha)
    
    # Compute the log ratio of the terms related to the proposal q
    grad_log_full_prop <- compute_grad_log_full_gamma_ibp(s_hat_prop, alpha_hat_prop, 
                                                            Gam, n, K, counts, 
                                                            a_s, b_s, a_alpha, b_alpha)
    
    log_ratio_q <- compute_log_ratio_q_gamma_ibp(s_hat_prop, alpha_hat_prop,
                                                   s_hat_curr, alpha_hat_curr, tau,
                                                   grad_log_full_curr, grad_log_full_prop)
    
    # Compute acceptance probability
    log_acc_prob <- log_ratio_full + log_ratio_q
    acc_prob <- min(1, exp(log_acc_prob))
    
    # Decide if accept or not the new parameter vector
    if (runif(1) < acc_prob){ # accept
      s <- exp(s_hat_prop)
      alpha <- exp(alpha_hat_prop)/(1 + exp(alpha_hat_prop))
    }
    
    
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


######################################################################


#' Metropolis-within-Gibbs sampler for BB + NegBin (Geometric prior on nstar)
#'
#' @param Z [integer] binary matrix of presence/absence (n x K- dimensional)
#' @param a_0 [numeric] initial value of a
#' @param b_0 [numeric] initial value of b
#' @param s_0 [numeric] initial value of s
#' @param alpha_0 [numeric] initial value of alpha 
#' @param p [numeric] hyperparameter of Geometric prior on a
#' @param r [numeric] 
#' @param t [numeric]
#' @param a_alpha [numeric]
#' @param b_alpha [numeric]
#' @param a_s [numeric]
#' @param b_s [numeric]
#' @param tau [numeric] MALA step-size
#' @param fixed [logical] vector indicating if the variable is fixed or has prior:
#' if fixed, it stays equal to initial value - order: a, b, alpha, s
#' @param S [integer] number of iterations for the MCMC algorithm
#' @param n_burnin [integer] number of iterations for the burn-in
#' @param thin [integer] thinning
#' @param seed [integer] seed
#'
#' @return
#' @import numDeriv
#' @import stats
#' @import MASS
#' @import matlib
#' @export
#'
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


######################################################################


#' Metropolis-within-Gibbs sampler for BB + NegBin (Geometric prior on nstar)
#'
#' @param Z [integer] binary matrix of presence/absence (n x K- dimensional)
#' @param a_0 [numeric] initial value of a
#' @param b_0 [numeric] initial value of b
#' @param s_0 [numeric] initial value of s
#' @param alpha_0 [numeric] initial value of alpha 
#' @param p [numeric] hyperparameter of Geometric prior on a
#' @param r [numeric] 
#' @param t [numeric]
#' @param a_alpha [numeric]
#' @param b_alpha [numeric]
#' @param a_s [numeric]
#' @param b_s [numeric]
#' @param tau [numeric] MALA step-size
#' @param fixed [logical] vector indicating if the variable is fixed or has prior:
#' if fixed, it stays equal to initial value - order: a, b, alpha, s
#' @param S [integer] number of iterations for the MCMC algorithm
#' @param n_burnin [integer] number of iterations for the burn-in
#' @param thin [integer] thinning
#' @param seed [integer] seed
#'
#' @return
#' @export
#'
#' @examples
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
    print(paste0("iteration: ", q))
    # pbar <- p_kmn_all_gamma_IBP(alpha, theta, M, n, b)
    # print(paste0("pbar: ", pbar))
    # kmn_chain[,q] <- rnbinom(M, a + Kn, pbar)
  }
  
  return (kmn_chain)
}

