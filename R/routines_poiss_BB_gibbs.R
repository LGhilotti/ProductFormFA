

#' Function the gradient of the log full-conditional for s_hat and alpha_bar_hat
#'
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
  
  dalpha_bar <- K + a_alpha + ea_bar*( (digamma(n+es+ea_bar) - digamma(es+ea_bar))*
                                         (- lambda* p1 - K) +
                                         sum(digamma(counts + ea_bar)) -
                                         K*digamma(1+ea_bar) - b_alpha)
  
  return (c(ds_hat, dalpha_bar))
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
compute_log_ratio_q <- function(s_hat_prop, alpha_bar_hat_prop,
                    s_hat_curr, alpha_bar_hat_curr, tau,
                    grad_log_full_curr, grad_log_full_prop){
  
  norm_prop <- norm(c(s_hat_curr, alpha_bar_hat_curr) - c(s_hat_prop, alpha_bar_hat_prop) -
    tau*grad_log_full_prop, type="2")^2
  
  norm_curr <- norm( c(s_hat_prop, alpha_bar_hat_prop) - c(s_hat_curr, alpha_bar_hat_curr) -
                      tau*grad_log_full_curr, type="2")^2
  
  res <- 1/(4*tau)*(norm_curr - norm_prop)
  
  return (res)
}


############################################################################


#' Metropolis-within-Gibbs sampler for BB + Poiss
#'
#' @param Z [integer] binary matrix of presence/absence (n x K- dimensional)
#' @param lambda_0 [numeric] initial value of lambda
#' @param alpha_bar_0 [numeric] initial value of alpha_bar 
#' @param s_0 [numeric] initial value of s
#' @param a_l [numeric] 
#' @param b_l [numeric]
#' @param a_alpha [numeric]
#' @param b_alpha [numeric]
#' @param a_s [numeric]
#' @param b_s [numeric]
#' @param tau [numeric] MALA step-size
#' @param S [integer] number of iterations for the MCMC algorithm
#' @param n_burnin [integer] number of iterations for the burn-in
#' @param thin [integer] thinning
#' @param seed [integer] seed
#'
#' @return
#' @export
#'
#' @examples
gibbs_sampler_poiss <- function(Z, 
                                lambda_0, alpha_bar_0, s_0,
                                a_l, b_l, a_alpha, b_alpha, a_s, b_s,
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
  
  ############## Gibbs-sampler ##########################
  
  # Define structure to store parameters along the iterations
  number_saved_iterations <- (S - n_burnin)/thin + 1
  lambda_vec <- vector(length = number_saved_iterations)
  alpha_bar_vec <- vector(length = number_saved_iterations)
  s_vec <- vector(length = number_saved_iterations)
  
  # Set initial values
  lambda <- lambda_0
  alpha_bar <- alpha_bar_0
  s <- s_0
  
  # index saved iterations (after burn-in and thinning satisfied)
  l <- 1
  
  for (q in 1:S){
    
    
    ################################################################
    ############# Draw lambda | Z, alpha_bar, s ###################
    ###############################################################
    
    t <- lgamma(s + n) + lgamma(s + alpha_bar) - lgamma(s) - lgamma(s+alpha_bar+n)
    
    lambda <- rgamma(1, shape = K + a_l ,
                     rate = b_l + 1 - exp(t))
    
    
    ################################################################
    ############# Draw (s, alpha_bar) | Z, lambda ##################
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
    
    s_hat_prop <- s_hat_curr + tau*grad_log_full[1] + sqrt(2*tau)*rnorm(1)
    alpha_bar_hat_prop <- alpha_bar_hat_curr + tau*grad_log_full[2] + sqrt(2*tau)*rnorm(1)
    
    ### Acceptance probability 
    # Compute the log ratio of the full-cond in prop point and curr point
    log_ratio_full <- compute_log_ratio_full_poiss(s_hat_prop, alpha_bar_hat_prop,
                                                   s_hat_curr, alpha_bar_hat_curr, lambda,
                                                   n, K, counts, a_s, b_s, a_alpha, b_alpha)
    
    # Compute the log ratio of the terms related to the proposal q
    grad_log_full_prop <- compute_grad_log_full_poiss(s_hat_prop, alpha_bar_hat_prop, lambda,
                                                      n, K, counts, a_s, b_s, a_alpha, b_alpha)
    
    log_ratio_q <- compute_log_ratio_q(s_hat_prop, alpha_bar_hat_prop,
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
      
      lambda_vec[l] <- lambda
      alpha_bar_vec[l] <- alpha_bar
      s_vec[l] <- s
      
      l <- l+1
    }
    
  }
  
  lambda_vec <- lambda_vec[1:(l-1)]
  alpha_bar_vec <- alpha_bar_vec[1:(l-1)]
  s_vec <- s_vec[1:(l-1)]
  
  return (list("lambda_vec" = lambda_vec, "alpha_bar_vec" = alpha_bar_vec, 
               "s_vec" = s_vec))
  
}




