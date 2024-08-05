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
      !all(class(initialization) == c("initialization", model)) ){
    stop("Prior and/or initialization not compatible")
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
  
  if (model == "GammaIBP_more_prior") {
    
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
                            sigq_alpha = sigq_alpha, sigq_s = sigq_s, S = S, n_burnin = n_burnin, thin = thin, seed = seed,
                            fix_a = FALSE, fix_b = FALSE)
    
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
  
  if (model == "GammaIBP_single_prior") {
    
    # Initialization of the chain
    alpha_0 <- initialization$alpha_0
    s_0 <- initialization$s_0
    # Hyperparameters
    a <- prior$a
    b <- prior$b
    a_alpha <- prior$a_alpha
    b_alpha <- prior$b_alpha
    a_s <- prior$a_s
    b_s <- prior$b_s
    # Additional MCMC parameters
    sigq_alpha <- mcmcparams$sigq_alpha
    sigq_s <- mcmcparams$sigq_s
    
    
    # Run the model
    res <- sampler_GammaIBP(Z = feature_matrix,
                            alpha_0 = alpha_0, s_0 = s_0, a_0 = a, b_0 = b,
                            a_alpha = a_alpha, b_alpha = b_alpha, a_s = a_s, b_s = b_s, 
                            q = 0.5, r = 1, t = 1,
                            sigq_alpha = sigq_alpha, sigq_s = sigq_s, S = S, n_burnin = n_burnin, thin = thin, seed = seed,
                            fix_a = TRUE, fix_b = TRUE)
    
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




#' Gibbs-type feature allocation models through Empirical Bayes (GibbsFA_eb): 
#' function to estimate the parameters via EB (maximizing the EFPF).
#'
#' @param feature_matrix A \code{n x K}-dimensional binary matrix of features
#' @param model Model to fit. Available models are \code{PoissonBB} (BB with Poisson(lambda) mixture),
#' \code{NegBinBB} (BB with NB(n0, mu0) mixture), \code{GammaIBP} (IBP with Gamma(a, b) mixture)
#' @param type Only "EFPF" 
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
  
}


