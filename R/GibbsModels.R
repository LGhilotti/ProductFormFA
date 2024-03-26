#' Gibbs-type feature allocation models (GibbsFA)
#'
#' @param feature_matrix A \code{n x K}-dimensional binary matrix of features
#' @param model Model to fit. Available models are \code{PoissonBB} (BB with Poisson(lambda) mixture),
#' \code{NegBinBB} (BB with NB(n0, mu0) mixture), \code{GammaIBP} (IBP with Gamma(a, b) mixture)
#' @param ... Additional parameters for hyperparameter elicitation and MCMC settins
#' 
#' @return An object of class \code{GibbsFA}
#'
#' @export
GibbsFA <- function(feature_matrix, model, hyperparams, initialization, mcmcparams, seed = 1234) {
  
  if (!all(class(hyperparams) == c("prior", model)) | 
      !all(class(initialization) == c("initialization", model)) |
      !all(class(mcmcparams) == c("mcmcparameters", model)) ){
    stop("Prior, initialization and/or MCMC parameters not compatible")
  }
  
  S <- mcmcparams$S
  n_burnin <- mcmcparams$n_burnin
  thin <- mcmcparams$thin
  
  if (model == "PoissonBB") {
    
    # Initialization of the chain
    alpha_bar_0 <- initialization$alpha_bar_0
    s_0 <- initialization$s_0
    # Hyperparameters
    a_alpha <- hyperparams$a_alpha
    b_alpha <- hyperparams$b_alpha
    a_s <- hyperparams$a_s
    b_s <- hyperparams$b_s
    lambda <- hyperparams$lambda
    # MCMC parameters
    tau <- mcmcparams$tau
    
    # Run the model
    res <- sampler_PoissonBB(Z = feature_matrix,
                             alpha_bar_0 = alpha_bar_0, s_0 = s_0,
                             a_alpha = a_alpha, b_alpha = b_alpha, a_s = a_s, b_s = b_s, lambda = lambda,
                             tau = tau, S = S, n_burnin = n_burnin, thin = thin, seed = seed)
    
    out <- list("feature_matrix" = feature_matrix,
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
    a_alpha <- hyperparams$a_alpha
    b_alpha <- hyperparams$b_alpha
    a_s <- hyperparams$a_s
    b_s <- hyperparams$b_s
    n0 <- hyperparams$n0
    mu0 <- hyperparams$mu0
    # MCMC parameters
    tau <- mcmcparams$tau
    
    # Run the model
    res <- sampler_NegBinBB(Z = feature_matrix,
                            alpha_bar_0 = alpha_bar_0, s_0 = s_0,
                            a_alpha = a_alpha, b_alpha = b_alpha, a_s = a_s, b_s = b_s, nstar = n0, p = mu0,
                            tau = tau, S = S, n_burnin = n_burnin, thin = thin, seed = seed)
    
    out <- list("feature_matrix" = feature_matrix,
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
    a_alpha <- hyperparams$a_alpha
    b_alpha <- hyperparams$b_alpha
    a_s <- hyperparams$a_s
    b_s <- hyperparams$b_s
    q <- hyperparams$q
    r <- hyperparams$r
    t <- hyperparams$t
    # MCMC parameters
    sigq_alpha <- mcmcparams$sigq_alpha
    sigq_s <- mcmcparams$sigq_s
    
    
    # Run the model
    res <- sampler_GammaIBP(Z = feature_matrix,
                            alpha_0 = alpha_0, s_0 = s_0, a_0 = a_0, b_0 = b_0,
                            a_alpha = a_alpha, b_alpha = b_alpha, a_s = a_s, b_s = b_s, q = q, r = r, t = t,
                            sigq_alpha = sigq_alpha, sigq_s = sigq_s, S = S, n_burnin = n_burnin, thin = thin, seed = seed)
    
    out <- list("feature_matrix" = feature_matrix,
                "alpha_chain" = res$alpha_chain, 
                "theta_chain" = res$s_chain - res$alpha_chain,
                "a_chain" = res$a_chain,
                "b_chain" = res$b_chain)
    
    class(out) <- c("GibbsFA", "GammaIBP")
    return(out)
  }
}