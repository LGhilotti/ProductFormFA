#' Gibbs-type feature allocation models (GibbsFA): function to fit the model (choose among 
#' PoissonBB, NegBinBB and GammaIBP).
#'
#' @param feature_matrix A \code{n x K}-dimensional binary matrix of features
#' @param model Model to fit. Available models are \code{PoissonBB} (BB with Poisson(lambda) mixture),
#' \code{NegBinBB} (BB with NB(n0, mu0) mixture), \code{GammaIBP} (IBP with Gamma(a, b) mixture)
#' @param prior Prior object for hyperparameter elicitation 
#' @param initialization Initialization object for parameters initialization
#' @param mcmcparams mcmcparameters object for MCMC setting
#' @param seed seed for fixing randomness
#' 
#' @return An object of class \code{GibbsFA, model}
#'
#' @export
#' 
GibbsFA <- function(feature_matrix, model, prior, initialization, mcmcparams, seed = 1234) {
  
  if (!all(class(prior) == c("prior", model)) | 
      !all(class(initialization) == c("initialization", model)) |
      !all(class(mcmcparams) == c("mcmcparameters", model)) ){
    stop("Prior, initialization and/or MCMC parameters not compatible")
  }
  
  # Remove NAs and 0s column from feature_matrix
  feature_matrix <- feature_matrix[, colSums(is.na(feature_matrix))==0]
  feature_matrix <- feature_matrix[, colSums(feature_matrix)!=0]
  
  # Extract MCMC setting of the chains
  S <- mcmcparams$S
  n_burnin <- mcmcparams$n_burnin
  thin <- mcmcparams$thin
  
  if (model == "PoissonBB") {
    
    # Initialization of the chain
    alpha_bar_0 <- initialization$alpha_bar_0
    s_0 <- initialization$s_0
    # Hyperparameters
    a_alpha <- prior$a_alpha
    b_alpha <- prior$b_alpha
    a_s <- prior$a_s
    b_s <- prior$b_s
    lambda <- prior$lambda
    # Additional MCMC parameters
    tau <- mcmcparams$tau
    
    # Run the model
    res <- sampler_PoissonBB(Z = feature_matrix,
                             alpha_bar_0 = alpha_bar_0, s_0 = s_0,
                             a_alpha = a_alpha, b_alpha = b_alpha, a_s = a_s, b_s = b_s, lambda = lambda,
                             tau = tau, S = S, n_burnin = n_burnin, thin = thin, seed = seed)
    
    out <- list("feature_matrix" = feature_matrix,
                "prior" = prior,
                "initialization" = initialization,
                "MCMCparameters" = mcmcparams,
                "alpha_chain" = - res$alpha_bar_chain, 
                "theta_chain" = res$alpha_bar_chain + res$s_chain)
    
    class(out) <- c("GibbsFA", "PoissonBB")
    return(out)
  }
  
  if (model == "NegBinBB") {
    
    # Initialization of the chain
    alpha_bar_0 <- initialization$alpha_bar_0
    s_0 <- initialization$s_0
    # Hyperparameters
    a_alpha <- prior$a_alpha
    b_alpha <- prior$b_alpha
    a_s <- prior$a_s
    b_s <- prior$b_s
    n0 <- prior$n0
    mu0 <- prior$mu0
    # Additional MCMC parameters
    tau <- mcmcparams$tau
    
    # Run the model
    res <- sampler_NegBinBB(Z = feature_matrix,
                            alpha_bar_0 = alpha_bar_0, s_0 = s_0,
                            a_alpha = a_alpha, b_alpha = b_alpha, a_s = a_s, b_s = b_s, nstar = n0, p = 1/(mu0/n0 + 1),
                            tau = tau, S = S, n_burnin = n_burnin, thin = thin, seed = seed)
    
    out <- list("feature_matrix" = feature_matrix,
                "prior" = prior,
                "initialization" = initialization,
                "MCMCparameters" = mcmcparams,
                "alpha_chain" = - res$alpha_bar_chain, 
                "theta_chain" = res$alpha_bar_chain + res$s_chain)
    
    class(out) <- c("GibbsFA", "NegBinBB")
    return(out)
  }
  
  if (model == "GammaIBP") {
    
    # Initialization of the chain
    alpha_0 <- initialization$alpha_0
    s_0 <- initialization$s_0
    a_0 <- initialization$a_0
    b_0 <- initialization$b_0
    # Hyperparameters
    a_alpha <- prior$a_alpha
    b_alpha <- prior$b_alpha
    a_s <- prior$a_s
    b_s <- prior$b_s
    q <- prior$q
    r <- prior$r
    t <- prior$t
    # Additional MCMC parameters
    sigq_alpha <- mcmcparams$sigq_alpha
    sigq_s <- mcmcparams$sigq_s
    
    
    # Run the model
    res <- sampler_GammaIBP(Z = feature_matrix,
                            alpha_0 = alpha_0, s_0 = s_0, a_0 = a_0, b_0 = b_0,
                            a_alpha = a_alpha, b_alpha = b_alpha, a_s = a_s, b_s = b_s, q = q, r = r, t = t,
                            sigq_alpha = sigq_alpha, sigq_s = sigq_s, S = S, n_burnin = n_burnin, thin = thin, seed = seed)
    
    out <- list("feature_matrix" = feature_matrix,
                "prior" = prior,
                "initialization" = initialization,
                "MCMCparameters" = mcmcparams,
                "alpha_chain" = res$alpha_chain, 
                "theta_chain" = res$s_chain - res$alpha_chain,
                "a_chain" = res$a_chain,
                "b_chain" = res$b_chain)
    
    class(out) <- c("GibbsFA", "GammaIBP")
    return(out)
  }
}



neg_log_EFPF_GibbsFA_R <- function(model, n, counts, par, known){
  
  par <- ifelse(is.na(known), par, known)
  
  if (model == "PoissonBB"){
    return (neg_log_EFPF_PoissonBB(n, counts, par))
  } 
  
  if (model == "NegBinBB"){
    # from parametrization (var_fct, mu0) need to find (n0,p)
    var_fct <- par[3]
    mu0 <- par[4]
    
    par[3] <- mu0/(var_fct -1) # this is n0
    par[4] <- 1/(mu0/par[3] + 1)
    return (neg_log_EFPF_NegBinBB(n, counts, par))
  }
  
  if (model == "NegBinBB_np"){
    return (neg_log_EFPF_NegBinBB(n, counts, par))
  }
  
  if (model == "GammaIBP"){
    return (neg_log_EFPF_GammaIBP(n, counts, par))
  }
  
}



neg_log_EFPF_BB_R <- function(n, counts, par, known){ # par: alpha, s, Nhat' = Nhat - k
  
  par <- ifelse(is.na(known), par, known)
  
  return(neg_log_EFPF_BB(n, counts, par))
  
}

neg_log_EFPF_IBP_R <- function(n, counts, par, known){ # par: alpha, s, Gamma
  
  par <- ifelse(is.na(known), par, known)
  
  return(neg_log_EFPF_IBP(n, counts, par))
  
}



mse_GammaIBP <- function(par, emp_knr, rmax, known){ # par = (alpha, s, a, b)
  
  par <- ifelse(is.na(known), par, known)
  
  alpha <- par[1]
  theta <- par[2] - par[1]
  a <- par[3]
  b <- par[4]

  return( sum((emp_knr[1:rmax] - ev_K_n_r_GammaIBP(alpha, theta, a, b, rmax))^2))
  
}

mse_PoissonBB <- function(par, emp_knr, rmax, known){ # par = (alpha, s, lambda)
  
  par <- ifelse(is.na(known), par, known)
  
  alpha <- par[1]
  theta <- par[2] - par[1]
  lambda <- par[3]

  return( sum((emp_knr[1:rmax] - ev_K_n_r_PoissonBB(alpha, theta, lambda, rmax))^2))
  
}


mse_NegBinBB <- function(par, emp_knr, rmax, known){ # par = (alpha, s, n0, mu0)
  
  par <- ifelse(is.na(known), par, known)
  
  alpha <- par[1]
  theta <- par[2] - par[1]
  n0 <- par[3]
  mu0 <- par[4]

  return( sum((emp_knr[1:rmax] - ev_K_n_r_NegBinBB(alpha, theta, n0, mu0, rmax))^2))
  
}



#' @export
mm_censored_fun <- function(x, n, emp_mean, emp_var){ # takes x = (a, b), where a = -alpha, b=alpha +theta
  
  a <- x[1]
  b <- x[2]
  beta_a_b <- beta(a,b)
  beta_a_b_n <- beta(a, b + n)
  
  ev_censored <- a/(beta_a_b - beta_a_b_n) * 
    (beta_a_b/(a + b) - beta_a_b_n/(a + b + n))
    
  y <- numeric(2)
  
  y[1] <- ev_censored - emp_mean
  y[2] <- beta_a_b/(beta_a_b - beta_a_b_n)*(a/(a+b))*((a+1)/(a+b+1)) +
    beta_a_b_n/(beta_a_b_n - beta_a_b)*(a/(a+b+n))*((a+1)/(a+b+n+1)) -
    ev_censored^2 - emp_var
  
  y
  
}

#' Gibbs-type feature allocation models through Empirical Bayes (GibbsFA_eb): 
#' function to estimate the parameters via EB, either by maximizing the EFPF
#' or method of moments (choose among PoissonBB, NegBinBB and GammaIBP).
#'
#' @param feature_matrix A \code{n x K}-dimensional binary matrix of features
#' @param model Model to fit. Available models are \code{PoissonBB} (BB with Poisson(lambda) mixture),
#' \code{NegBinBB} (BB with NB(n0, mu0) mixture), \code{GammaIBP} (IBP with Gamma(a, b) mixture)
#' @param type Either "EFPF" (for all models),"MM_biased", "MM_censored" (only for mixtures of BBs),
#' "MM" (only for GammaIBP)
#' @param seed seed for fixing randomness
#' 
#' @return An object of class \code{GibbsFA, model_eb}
#'
#' @import nleqslv
#' @export
#' 
GibbsFA_eb <- function(feature_matrix, model, type, seed = 1234, 
                       eb_params = NULL, 
                       Nhat_MM = NULL, var_fct = NULL, 
                       var_GammaIBP = NULL, ...) {
  
  if (type == "EFPF"){
    
    # Remove NAs and 0s column from feature_matrix
    feature_matrix <- feature_matrix[, colSums(is.na(feature_matrix))==0]
    feature_matrix <- feature_matrix[, colSums(feature_matrix)!=0]
    counts <- colSums(feature_matrix)
    n <- nrow(feature_matrix)
    K <- ncol(feature_matrix)
    
    if (all(class(eb_params) == c("eb_params", "BB"))){ # optimizing alpha, theta, Nhat_prime
  
      if (model == "PoissonBB" | model == "NegBinBB") {
        
        eb_init <- eb_params$init # this contains Nhat_prime, NOT Nhat
        eb_known <- eb_params$known
        
        res <- nlminb(
          start = eb_init, objective =  neg_log_EFPF_BB_R,  
          n = n, counts = counts, known = eb_known, lower = c(-Inf, 1e-5, 1e-5), upper = c(-1e-5, Inf, Inf)
        )
        
        Nhat_res <- round(unname(res$par["Nhat_prime"]) + K)
        
        if (model == "PoissonBB"){
          out <- list("feature_matrix" = feature_matrix,
                      "eb_params" = eb_params,
                      "alpha" = unname(res$par["alpha"]), 
                      "theta" = unname(res$par["s"] - res$par["alpha"]),
                      "lambda" = Nhat_res,
                      "fun_value" = res$objective
          )
          
          class(out) <- c("GibbsFA", "PoissonBB_eb")
          return(out)
        }
        
        if (model == "NegBinBB"){
          out <- list("feature_matrix" = feature_matrix,
                      "eb_params" = eb_params,
                      "alpha" = unname(res$par["alpha"]), 
                      "theta" = unname(res$par["s"] - res$par["alpha"]),
                      "var_fct" = var_fct,
                      "n0" = Nhat_res/(var_fct - 1),
                      "mu0" = Nhat_res,
                      "fun_value" = res$objective
          )
         
          class(out) <- c("GibbsFA", "NegBinBB_eb")
          return(out)
        }
        
      }
    } 
    
    if (all(class(eb_params) == c("eb_params", "IBP"))){ # optimizing alpha, theta, Gamma
      
      if (model == "GammaIBP"){ # we always optimize all the parameters, never fix them
        
        eb_init <- eb_params$init 
        eb_known <- eb_params$known
        
        res <- nlminb(
          start = eb_init, objective =  neg_log_EFPF_IBP_R,  
          n = n, counts = counts, known = eb_known, 
          lower = c(1e-5, 1e-5, 1e-5), upper = c(1 - 1e-5, Inf, Inf)
        )
        
        Gamma_prior_mean <- unname(res$par["Gamma"])
        
        out <- list("feature_matrix" = feature_matrix,
                    "eb_params" = eb_params,
                    "alpha" = unname(res$par["alpha"]), 
                    "theta" = unname(res$par["s"] - res$par["alpha"]),
                    "var" = var_GammaIBP,
                    "a" = Gamma_prior_mean^2 / var_GammaIBP,
                    "b" = Gamma_prior_mean / var_GammaIBP,
                    "fun_value" = res$objective
        )
        
        class(out) <- c("GibbsFA", "GammaIBP_eb")
        return(out)
        
      }
      
    }
      
      
    if (!all(class(eb_params) == c("eb_params", model))  ){
      stop("Starting point/known parameters for optimization not compatible")
    }
    
    if (model == "PoissonBB") {
      
      # Initialization of the optimization
      eb_init <- eb_params$init
      eb_known <- eb_params$known
      
      res <- nlminb(
        start = eb_init, objective =  neg_log_EFPF_GibbsFA_R, model = "PoissonBB", 
        n = n, counts = counts, known = eb_known, lower = c(-Inf, 1e-5, 1e-5), upper = c(-1e-5, Inf, Inf)
      )
      
      
      out <- list("feature_matrix" = feature_matrix,
                  "eb_params" = eb_params,
                  "alpha" = unname(res$par["alpha"]), 
                  "theta" = unname(res$par["s"] - res$par["alpha"]),
                  "lambda" = eb_known[["lambda"]],
                  "fun_value" = res$objective
      )
      
      class(out) <- c("GibbsFA", "PoissonBB_eb")
      return(out)
    }
    
    if (model == "NegBinBB") {
      
      # Initialization of the optimization
      eb_init <- eb_params$init
      eb_known <- eb_params$known
      
      res <- nlminb(
        start = eb_init, objective = neg_log_EFPF_GibbsFA_R, model = "NegBinBB", 
        n = n, counts = counts, known = eb_known,
        lower = c(-Inf, 1e-5, 1 + 1e-5, 1e-5), upper = c(-1e-5, Inf, Inf, Inf)
      )
      
      
      out <- list("feature_matrix" = feature_matrix,
                  "eb_params" = eb_params,
                  "alpha" = unname(res$par["alpha"]), 
                  "theta" = unname(res$par["s"] - res$par["alpha"]),
                  "var_fct" = eb_known[["var_fct"]],
                  "n0" = eb_known[["mu0"]]/(eb_known[["var_fct"]] - 1),
                  "mu0" = eb_known[["mu0"]],
                  "fun_value" = res$objective
      )
      
      class(out) <- c("GibbsFA", "NegBinBB_eb")
      return(out)
    }
    
    if (model == "GammaIBP") {
      stop("not implemented")
    }
    
  }
  
  
  # if (type == "EFPF_wrong"){
  #   
  #   if (!all(class(eb_params) == c("eb_params", model))  ){
  #     stop("Starting point/known parameters for optimization not compatible")
  #   }
  #   
  #   # Remove NAs and 0s column from feature_matrix
  #   feature_matrix <- feature_matrix[, colSums(is.na(feature_matrix))==0]
  #   feature_matrix <- feature_matrix[, colSums(feature_matrix)!=0]
  #   counts <- colSums(feature_matrix)
  #   n <- nrow(feature_matrix)
  #   
  #   if (model == "PoissonBB") {
  #     
  #     # Initialization of the optimization
  #     eb_init <- eb_params$init
  #     eb_known <- eb_params$known
  #     
  #     res <- nlminb(
  #       start = eb_init, objective =  neg_log_EFPF_GibbsFA_R, model = "PoissonBB", 
  #       n = n, counts = counts, known = eb_known, lower = c(-Inf, 1e-5, 1e-5), upper = c(-1e-5, Inf, Inf)
  #     )
  #     
  #     
  #     out <- list("feature_matrix" = feature_matrix,
  #                 "eb_params" = eb_params,
  #                 "alpha" = unname(res$par["alpha"]), 
  #                 "theta" = unname(res$par["s"] - res$par["alpha"]),
  #                 "lambda" = unname(res$par["lambda"]),
  #                 "fun_value" = res$objective
  #     )
  #     
  #     class(out) <- c("GibbsFA", "PoissonBB_eb")
  #     return(out)
  #   }
  #   
  #   if (model == "NegBinBB") {
  #     
  #     # Initialization of the optimization
  #     eb_init <- eb_params$init
  #     eb_known <- eb_params$known
  #     
  #     res <- nlminb(
  #       start = eb_init, objective = neg_log_EFPF_GibbsFA_R, model = "NegBinBB", 
  #       n = n, counts = counts, known = eb_known,
  #       lower = c(-Inf, 1e-5, 1 + 1e-5, 1e-5), upper = c(-1e-5, Inf, Inf, Inf)
  #     )
  #     
  #     
  #     out <- list("feature_matrix" = feature_matrix,
  #                 "eb_params" = eb_params,
  #                 "alpha" = unname(res$par["alpha"]), 
  #                 "theta" = unname(res$par["s"] - res$par["alpha"]),
  #                 "var_fct" = unname(res$par["var_fct"]),
  #                 "n0" = unname(res$par["mu0"]/(res$par["var_fct"] - 1)),
  #                 "mu0" = unname(res$par["mu0"]),
  #                 "fun_value" = res$objective
  #     )
  #     
  #     class(out) <- c("GibbsFA", "NegBinBB_eb")
  #     return(out)
  #   }
  #   
  #   
  #   if (model == "NegBinBB_np") {
  #     
  #     # Initialization of the optimization
  #     eb_init <- eb_params$init
  #     eb_known <- eb_params$known
  #     
  #     res <- optim(
  #       par = eb_init, fn = neg_log_EFPF_GibbsFA_R, model = "NegBinBB_np", 
  #       n = n, counts = counts, known = eb_known,
  #       method = "L-BFGS-B", lower = c(-Inf, 1e-5, 1e-5, 1e-5), upper = c(-1e-5, Inf, Inf, 1 - 1e-5)
  #     )
  #     
  #     
  #     out <- list("feature_matrix" = feature_matrix,
  #                 "eb_params" = eb_params,
  #                 "alpha" = unname(res$par["alpha"]), 
  #                 "theta" = unname(res$par["s"] - res$par["alpha"]),
  #                 "n0" = unname(res$par["n0"]),
  #                 "p" = unname(res$par["p"])
  #     )
  #     
  #     class(out) <- c("GibbsFA", "NegBinBB_np_eb")
  #     return(out)
  #   }
  #   
  #   
  #   if (model == "GammaIBP") {
  #     
  #     # Initialization of the optimization
  #     eb_init <- eb_params$init
  #     eb_known <- eb_params$known
  #     
  #     res <- optim(
  #       par = eb_init, fn = neg_log_EFPF_GibbsFA_R, model = "GammaIBP", 
  #       n = n, counts = counts, known = eb_known,
  #       method = "L-BFGS-B", lower = c(1e-5, 1e-5, 1e-5, 1e-5), upper = c(1 - 1e-5, Inf, Inf, Inf)
  #     )
  #     
  #     out <- list("feature_matrix" = feature_matrix,
  #                 "eb_params" = eb_params,
  #                 "alpha" = unname(res$par["alpha"]), 
  #                 "theta" = unname(res$par["s"] - res$par["alpha"]),
  #                 "a" = unname(res$par["a"]),
  #                 "b" = unname(res$par["b"])
  #     )
  #     
  #     class(out) <- c("GibbsFA", "GammaIBP_eb")
  #     return(out)
  #   }
  # }
  # 
  # 
  # 
  # if (type == "MM_biased"){
  #   
  #   if (model == "GammaIBP"){
  #     stop("GammaIBP has not the MM_biased based parameter elicitation")
  #   }
  #   
  #   # Remove NAs and 0s column from feature_matrix
  #   feature_matrix <- feature_matrix[, colSums(is.na(feature_matrix))==0]
  #   feature_matrix <- feature_matrix[, colSums(feature_matrix)!=0]
  #   counts <- colSums(feature_matrix)
  #   n <- nrow(feature_matrix)
  #   K <- ncol(feature_matrix)
  #   
  #   # Compute empirical mean and variance of observed pis
  #   emp_pis <- counts/n
  #   emp_mean <- mean(emp_pis)
  #   emp_var <- var(emp_pis)
  #   
  #   # Solve the system given by the moments equations
  #   theta_MM <- emp_mean*(1 - emp_mean)/emp_var - 1
  #   alpha_MM <- - emp_mean*theta_MM
  #   if (is.null(Nhat_MM)){
  #     Nhat_MM <- K/(feature_fraction(n, alpha_MM, theta_MM))
  #   }
  #   
  #   out <- list("feature_matrix" = feature_matrix,
  #               "alpha" = alpha_MM, 
  #               "theta" = theta_MM)
  #   
  #   if (model == "PoissonBB") {
  #   
  #     out[["lambda"]] <- Nhat_MM
  #     
  #     class(out) <- c("GibbsFA", "PoissonBB_eb")
  #     return(out)
  #   }
  #   
  #   if (model == "NegBinBB") {
  #     
  #     out[["mu0"]] <- Nhat_MM
  #     out[["n0"]] <- Nhat_MM/(var_fct - 1)
  #       
  #     class(out) <- c("GibbsFA", "NegBinBB_eb")
  #     return(out)
  #   }
  #   
  # }
  #   
  # 
  # if (type == "MM_censored"){
  #   
  #   if (model == "GammaIBP"){
  #     stop("GammaIBP has not the MM_censored based parameter elicitation")
  #   }
  #   
  #   # Remove NAs and 0s column from feature_matrix
  #   feature_matrix <- feature_matrix[, colSums(is.na(feature_matrix))==0]
  #   feature_matrix <- feature_matrix[, colSums(feature_matrix)!=0]
  #   counts <- colSums(feature_matrix)
  #   n <- nrow(feature_matrix)
  #   K <- ncol(feature_matrix)
  #   
  #   # MM estimates are just used as initialization of the system solver algorithm
  #   emp_pis <- counts/n
  #   emp_mean <- mean(emp_pis)
  #   emp_var <- var(emp_pis)
  #   theta_MM <- emp_mean*(1 - emp_mean)/emp_var - 1
  #   alpha_MM <- - emp_mean*theta_MM
  #   
  #   # Solve the non-linear system of equations for alpha and theta 
  #   x_start <- c(- alpha_MM, alpha_MM + theta_MM)
  #   sol <- nleqslv(x_start, mm_censored_fun, control=list(btol=.01),
  #                  n = n, emp_mean = emp_mean, emp_var = emp_var)
  #   
  #   alpha_MM_censored <- - sol$x[1]
  #   theta_MM_censored <- sol$x[1] + sol$x[2]
  #   
  #   if (is.null(Nhat_MM)){
  #     Nhat_MM <- K/(feature_fraction(n, alpha_MM_censored, theta_MM_censored))
  #   }    
  #   
  #   out <- list("feature_matrix" = feature_matrix,
  #               "alpha" = alpha_MM_censored, 
  #               "theta" = theta_MM_censored)
  #   
  #   if (model == "PoissonBB") {
  #     
  #     out[["lambda"]] <- Nhat_MM
  #     
  #     class(out) <- c("GibbsFA", "PoissonBB_eb")
  #     return(out)
  #   }
  #   
  #   if (model == "NegBinBB") {
  #     
  #     out[["mu0"]] <- Nhat_MM
  #     out[["n0"]] <- Nhat_MM/(var_fct - 1)
  #     
  #     class(out) <- c("GibbsFA", "NegBinBB_eb")
  #     return(out)
  #   }
  #   
  # }
  # 
  # if (type == "MM_knr"){
  #   
  #   if (!all(class(eb_params) == c("eb_params", model))  ){
  #     stop("Starting point/known parameters for optimization not compatible")
  #   }
  #   
  #   # Remove NAs and 0s column from feature_matrix
  #   feature_matrix <- feature_matrix[, colSums(is.na(feature_matrix))==0]
  #   feature_matrix <- feature_matrix[, colSums(feature_matrix)!=0]
  #   n <- nrow(feature_matrix)
  #   K <- ncol(feature_matrix)
  #   
  #   if (rmax == 0){
  #     rmax <- n
  #   }
  #   
  #   if (rmax > n){
  #     stop("rmax has to be less or equal to n")
  #   }
  #     
  #   emp_knr <- K_n_r(feature_matrix)[[paste0("N = ", n)]]
  #   
  #   if (model == "GammaIBP"){
  #     
  #     # Initialization of the optimization
  #     eb_init <- eb_params$init
  #     eb_known <- eb_params$known
  #     
  #     res <- nlminb(
  #       start = eb_init, objective = mse_GammaIBP, 
  #       emp_knr = emp_knr, rmax = rmax, known = eb_known,
  #       lower = c(1e-5, 1e-5, 1e-5, 1e-5), upper = c(1 - 1e-5, Inf, Inf, Inf)
  #     )
  #     
  #     out <- list("feature_matrix" = feature_matrix,
  #                 "eb_params" = eb_params,
  #                 "alpha" = unname(res$par["alpha"]), 
  #                 "theta" = unname(res$par["s"] - res$par["alpha"]),
  #                 "a" = unname(res$par["a"]),
  #                 "b" = unname(res$par["b"]),
  #                 "rmax" = rmax,
  #                 "fun_value" = res$objective
  #     )
  #     
  #     class(out) <- c("GibbsFA", "GammaIBP_eb")
  #     return(out)
  #     
  #   }
  #   
  #   if (model == "PoissonBB"){
  #     
  #     # Initialization of the optimization
  #     eb_init <- eb_params$init
  #     eb_known <- eb_params$known
  #     
  #     res <- nlminb(
  #       start = eb_init, objective = mse_PoissonBB, 
  #       emp_knr = emp_knr, rmax = rmax, known = eb_known,
  #       lower = c(-Inf, 1e-5, 1e-5), upper = c(-1e-5, Inf, Inf)
  #     )
  #     
  #     out <- list("feature_matrix" = feature_matrix,
  #                 "eb_params" = eb_params,
  #                 "alpha" = unname(res$par["alpha"]), 
  #                 "theta" = unname(res$par["s"] - res$par["alpha"]),
  #                 "lambda" = unname(res$par["lambda"]),
  #                 "rmax" = rmax,
  #                 "fun_value" = res$objective
  #     )
  #     
  #     class(out) <- c("GibbsFA", "PoissonBB_eb")
  #     return(out)
  #     
  #   }
  #   
  #   if (model == "NegBinBB"){
  #     
  #     # Initialization of the optimization
  #     eb_init <- eb_params$init
  #     eb_known <- eb_params$known
  #     
  #     res <- nlminb(
  #       start = eb_init, objective = mse_NegBinBB, 
  #       emp_knr = emp_knr, rmax = rmax, known = eb_known,
  #       lower = c(-Inf, 1e-5, 1 + 1e-5, 1e-5), upper = c(-1e-5, Inf, Inf, Inf)
  #     )
  #     
  #     out <- list("feature_matrix" = feature_matrix,
  #                 "eb_params" = eb_params,
  #                 "alpha" = unname(res$par["alpha"]), 
  #                 "theta" = unname(res$par["s"] - res$par["alpha"]),
  #                 "var_fct" = unname(res$par["var_fct"]),
  #                 "n0" = unname(res$par["mu0"]/(res$par["var_fct"] - 1)),
  #                 "mu0" = unname(res$par["mu0"]),
  #                 "rmax" = rmax,
  #                 "fun_value" = res$objective
  #     )
  #     
  #     class(out) <- c("GibbsFA", "NegBinBB_eb")
  #     return(out)
  #     
  #   }
  #   
  # }
}



#' Plot for both the mixtures of Beta-Bernoulli and the Gamma IBP
#'
#' @param x An object of class \code{GibbsFA}.
#' @param type Type of plot. Available options are:  \code{"rarefaction"}, \code{"extrapolation"} and \code{"richness"}.
#' Type \code{"richness"} is available only for mixtures of Beta-Bernoulli.
#' @param bw bandwidth for sample-based plot of richness.
#' @param n_reorderings number of reorderings for accumularion curve.
#' @param M additional sample to predict. Valid only for \code{type = "extrapolation"}. By default, it is equal to the sample size.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details The method produces summary plots for the GibbsFA output
#'
#' @import ggplot2
#' @import scales
#' @import ggpubr
#' @import dplyr
#' @export
#'
plot.GibbsFA <- function(x, type, bw = 1, n_reorderings = 1, M = NULL, seed = 1234, ...){
  
  if (!type %in% c("richness", "rarefaction", "extrapolation")){
    stop("Error: invalid type of plot")
  }
  
  n <- nrow(x$feature_matrix)
  Kn <- ncol(x$feature_matrix)
  
  
  if (type == "richness"){
    
    richness <- total_richness(object = x)
    
    if (class(x)[2] == "PoissonBB_eb"){
      
      lambda <- richness$lambda_post
      lb <- qpois(0.025, lambda, lower.tail = TRUE, log.p = FALSE)
      ub <- qpois(0.975, lambda, lower.tail = TRUE, log.p = FALSE)
      y <- dpois(lb:ub, lambda = lambda)
      
      p <- ggplot(data.frame(x= Kn + (lb:ub), y = y), aes(x,y)) +
        geom_point() +
        theme_light() +
        theme(legend.position = "top") +
        scale_y_continuous(breaks = pretty_breaks()) +
        xlab("# distinct features") + rremove("ylab") +
        facet_wrap(~"Richness") +
        theme(aspect.ratio = 1)
      
      
    } else if (class(x)[2] == "NegBinBB_eb"){

      n0_post <- richness$n0_post
      mu0_post <- richness$mu0_post
      p_post <- 1/(mu0_post/n0_post + 1)
      lb <- qnbinom(0.025, size = n0_post, prob = p_post, lower.tail = TRUE, log.p = FALSE)
      ub <- qnbinom(0.975, size = n0_post, prob = p_post, lower.tail = TRUE, log.p = FALSE)
      y <- dnbinom(lb:ub, size = n0_post, prob = p_post)
        
      p <- ggplot(data.frame(x= Kn + (lb:ub), y = y), aes(x,y)) +
        geom_point() +
        theme_light() +
        theme(legend.position = "top") +
        scale_y_continuous(breaks = pretty_breaks()) +
        xlab("# distinct features") + rremove("ylab") +
        facet_wrap(~"Richness") +
        theme(aspect.ratio = 1)
    
    } else {
    
      df <- data.frame("richness" = richness, "model" = class(x)[2])
      lb <- quantile(richness, prob = 0.025)
      ub <- quantile(richness, prob = 0.975)
      
      p <- ggplot(df) +
        stat_density(aes(x=richness), geom="line",position="identity", bw = bw) +
        theme_light() +
        theme(legend.position = "top") +
        scale_y_continuous(breaks = pretty_breaks()) +
        xlab("# distinct features") + rremove("ylab") +
        scale_x_continuous(limits = c(lb-1, ub+1)) +
        facet_wrap(~"Richness") +
        theme(aspect.ratio = 1)
    }
      
  } else if (type == "rarefaction"){
    
    rare_list <- rarefaction(object = x, seed = seed)

    # Extract accumulation curve of the observed sample (or average accumulation)
    accum <- rarefaction(x$feature_matrix, n_reorderings = n_reorderings, seed = seed)
    
    if (class(x)[2] == "PoissonBB_eb"){
      
      lambda_post <- rare_list$lambda_post
      lb_bands <- qpois(0.025, lambda_post, lower.tail = TRUE, log.p = FALSE)
      ub_bands <- qpois(0.975, lambda_post, lower.tail = TRUE, log.p = FALSE)
      means <- lambda_post
      df <- data.frame("means" = c(0,means), 
                       "lb_bands" = c(0,lb_bands),
                       "ub_bands" = c(0,ub_bands),
                       "accum" = c(0,accum))
      
      p <- ggplot(df, aes(x = 0:n, y = means)) +
        geom_line(linetype = "solid", color = "red" , linewidth = 0.9) +
        #geom_ribbon(aes(ymin = lb_bands, ymax = ub_bands), color = "red" , linewidth = 0.8, alpha = 0.1) +
        geom_point( aes(x = 0:n, y = accum), color="black", shape = 1) +
        xlab("# observations") + ylab("# distinct features") + 
        theme_light() + 
        theme(legend.position = "top") +
        scale_y_continuous(breaks = pretty_breaks()) +
        scale_x_continuous(breaks = pretty_breaks()) +
        facet_wrap(~"Rarefaction and accumulation") +
        theme(aspect.ratio = 1)
      
      
    } else if (class(x)[2] %in% c("NegBinBB_eb", "GammaIBP_eb")){
      
      n0_post <- rare_list$n0_post
      mu0_post <- rare_list$mu0_post
      p_post <- 1/(mu0_post/n0_post + 1)
      lb_bands <- qnbinom(0.025, size = n0_post, prob = p_post, lower.tail = TRUE, log.p = FALSE)
      ub_bands <- qnbinom(0.975, size = n0_post, prob = p_post, lower.tail = TRUE, log.p = FALSE)
      means <- mu0_post
      df <- data.frame("means" = c(0,means), 
                       "lb_bands" = c(0,lb_bands),
                       "ub_bands" = c(0,ub_bands),
                       "accum" = c(0,accum))
      
      p <- ggplot(df, aes(x = 0:n, y = means)) +
        geom_line(linetype = "solid", color = "red" , linewidth = 0.9) +
        #geom_ribbon(aes(ymin = lb_bands, ymax = ub_bands), color = "red" , linewidth = 0.8, alpha = 0.1) +
        geom_point( aes(x = 0:n, y = accum), color="black", shape = 1) +
        xlab("# observations") + ylab("# distinct features") + 
        theme_light() + 
        theme(legend.position = "top") +
        scale_y_continuous(breaks = pretty_breaks()) +
        scale_x_continuous(breaks = pretty_breaks()) +
        facet_wrap(~"Rarefaction and accumulation") +
        theme(aspect.ratio = 1)
      
    } else { # Prior approach
      
      ci_mat <- as.matrix(bind_rows(lapply(rare_list, quantile, probs = c(0.025,0.975) )))
      mean_vec <- as.vector(unlist(lapply(rare_list, mean)))
      
      df <- data.frame("means" = c(0,mean_vec),
                       "lbs" = c(0,ci_mat[,1]),
                       "ubs" = c(0,ci_mat[,2]),
                       "accum" = c(0, accum))
      
      p <- ggplot(df, aes(x = 0:n, y = means)) +
        geom_line(linetype = "solid", color = "red" , linewidth = 0.9) +
        #geom_ribbon(aes(ymin = lbs, ymax = ubs), color = "red" , linewidth = 0.8, alpha = 0.1) +
        geom_point( aes(x = 0:n, y = accum), color="black", shape = 1) +
        xlab("# observations") + ylab("# distinct features") + 
        theme_light() + 
        theme(legend.position = "top") +
        scale_y_continuous(breaks = pretty_breaks()) +
        scale_x_continuous(breaks = pretty_breaks()) +
        facet_wrap(~"Rarefaction and accumulation") +
        theme(aspect.ratio = 1)
      
    }
    
  
  } else if (type == "extrapolation"){
    
    if(is.null(M)){
      M = n
    }
    
    extr_list <- extrapolation(object = x, M = M, seed = seed)
    
    # Extract accumulation curve of the observed sample (or average accumulation)
    accum <- rarefaction(x$feature_matrix, n_reorderings = n_reorderings, seed = seed)
    
    if (class(x)[2] == "PoissonBB_eb"){
      
      lambda_post <- extr_list$lambda_post
      lb_bands <- qpois(0.025, lambda_post, lower.tail = TRUE, log.p = FALSE)
      ub_bands <- qpois(0.975, lambda_post, lower.tail = TRUE, log.p = FALSE)
      means <- lambda_post
      df <- data.frame("means" = c(Kn,Kn + means), 
                       "lb_bands" = c(Kn,Kn + lb_bands),
                       "ub_bands" = c(Kn,Kn + ub_bands))
      
      p <- ggplot(df, aes(x = n:(n+M), y = means)) +
        geom_line(linetype = "dashed", color = "red" , linewidth = 0.9) +
        geom_ribbon(aes(ymin = lb_bands, ymax = ub_bands), color = "red" , linewidth = 0.8, alpha = 0.1) +
        geom_point( data = data.frame(x = 0:n, y = c(0,accum)), aes(x = x, y = y),
                    color="black", shape = 1) +
        geom_vline(aes(xintercept = n), linetype="dashed", color = "grey") +
        xlab("# observations") + ylab("# distinct features") + 
        theme_light() + 
        theme(legend.position = "top") +
        scale_y_continuous(breaks = pretty_breaks()) +
        scale_x_continuous(breaks = pretty_breaks()) +
        facet_wrap(~"Extrapolation and accumulation") +
        theme(aspect.ratio = 1)
      
      
    } else if (class(x)[2] %in% c("NegBinBB_eb", "GammaIBP_eb")){
      
      n0_post <- extr_list$n0_post
      mu0_post <- extr_list$mu0_post
      p_post <- 1/(mu0_post/n0_post + 1)
      lb_bands <- qnbinom(0.025, size = n0_post, prob = p_post, lower.tail = TRUE, log.p = FALSE)
      ub_bands <- qnbinom(0.975, size = n0_post, prob = p_post, lower.tail = TRUE, log.p = FALSE)
      means <- mu0_post
      df <- data.frame("means" = c(Kn,Kn + means), 
                       "lb_bands" = c(Kn,Kn + lb_bands),
                       "ub_bands" = c(Kn,Kn + ub_bands))
      
      p <- ggplot(df, aes(x = n:(n+M), y = means)) +
        geom_line(linetype = "dashed", color = "red" , linewidth = 0.9) +
        geom_ribbon(aes(ymin = lb_bands, ymax = ub_bands), color = "red" , linewidth = 0.8, alpha = 0.1) +
        geom_point( data = data.frame(x = 0:n, y = c(0,accum)), aes(x = x, y = y),
                    color="black", shape = 1) +
        geom_vline(aes(xintercept = n), linetype="dashed", color = "grey") +
        xlab("# observations") + ylab("# distinct features") + 
        theme_light() + 
        theme(legend.position = "top") +
        scale_y_continuous(breaks = pretty_breaks()) +
        scale_x_continuous(breaks = pretty_breaks()) +
        facet_wrap(~"Extrapolation and accumulation") +
        theme(aspect.ratio = 1)
      
    } else { # Prior approach
      
      ci_mat <- as.matrix(bind_rows(lapply(extr_list, quantile, probs = c(0.025,0.975) )))
      mean_vec <- as.vector(unlist(lapply(extr_list, mean)))
      
      df <- data.frame("means" = c(Kn,mean_vec),
                       "lbs" = c(Kn,ci_mat[,1]),
                       "ubs" = c(Kn,ci_mat[,2]))
      
      p <- ggplot(df, aes(x = n:(n+M), y = means)) +
        geom_line(linetype = "dashed", color = "red" , linewidth = 0.9) +
        geom_ribbon(aes(ymin = lbs, ymax = ubs), color = "red" , linewidth = 0.8, alpha = 0.1) +
        geom_point( data = data.frame(x = 0:n, y = c(0,accum)), aes(x = x, y = y), 
                    color="black", shape = 1) +
        geom_vline(aes(xintercept = n), linetype="dashed", color = "grey") +
        xlab("# observations") + ylab("# distinct features") + 
        theme_light() + 
        theme(legend.position = "top") +
        scale_y_continuous(breaks = pretty_breaks()) +
        scale_x_continuous(breaks = pretty_breaks()) +
        facet_wrap(~"Extrapolation and accumulation") +
        theme(aspect.ratio = 1)
    }
    
  
    
  }
  
  return(p)
  
}