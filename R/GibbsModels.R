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
                            a_alpha = a_alpha, b_alpha = b_alpha, a_s = a_s, b_s = b_s, nstar = n0, p = mu0,
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




#' Plot for both the mixtures of Beta-Bernoulli and the Gamma IBP
#'
#' @param x An object of class \code{GibbsFA}.
#' @param M additional sample to predict. Valid only for \code{type = "extrapolation"}. By default, it is equal to the sample size.
#' @param type Type of plot. Available options are:  \code{"rarefaction"}, \code{"extrapolation"} and \code{"richness"}.
#' Type \code{"richness"} is available only for mixtures of Beta-Bernoulli.
#' @param plot_sample Whether to add the sample-based rarefaction curve. Available for \code{type = "rarefaction"}.
#' @param ... Further arguments passed to or from other methods.
#'
#' @details The method produces summary plots for the GibbsFA output
#'
#' @export
#'
plot.GibbsFA <- function(x, type, plot_sample = TRUE, M = NULL, seed = 1234, ...){
  
  if (!type %in% c("richness", "rarefaction", "extrapolation")){
    stop("Error: invalid type of plot")
  }
  
  n <- nrow(x$feature_matrix)
  Kn <- ncol(x$feature_matrix)
  
  
  if (type == "richness"){
    
    richness <- total_richness(object = x)
    
    df <- data.frame("richness" = richness, "model" = class(x)[2])
    
    p <- ggplot(df) +
      stat_density(aes(x=richness), geom="line",position="identity") +
      theme_light() +
      theme(legend.position = "top") +
      scale_y_continuous(breaks = pretty_breaks()) +
      xlab("# distinct features") + rremove("ylab") +
      facet_wrap(~"Richness") +
      theme(aspect.ratio = 1)
    
  } else if (type == "rarefaction"){
    
    rare_list <- rarefaction(object = x, seed = seed)
    
    ci_mat <- as.matrix(bind_rows(lapply(rare_list, quantile, probs = c(0.025,0.975) )))
    mean_vec <- as.vector(unlist(lapply(rare_list, mean)))
    
    if (plot_sample == TRUE){
      
      # Extract accumulation curve of the observed sample (or average accumulation)
      accum <- rarefaction(x$feature_matrix)
      
      df <- data.frame("means" = c(0,mean_vec),
                       "lbs" = c(0,ci_mat[,1]),
                       "ubs" = c(0,ci_mat[,2]),
                       "accum" = c(0, accum))
      
      p <- ggplot(df, aes(x = 0:n, y = means)) +
        geom_line(linetype = "dashed", color = "red" , linewidth = 0.9) +
        geom_ribbon(aes(ymin = lbs, ymax = ubs), color = "red" , linewidth = 0.8, alpha = 0.1) +
        geom_point( aes(x = 0:n, y = accum), color="black", shape = 1) +
        xlab("# observations") + ylab("# distinct features") + 
        theme_light() + 
        theme(legend.position = "top") +
        scale_y_continuous(breaks = pretty_breaks()) +
        scale_x_continuous(breaks = pretty_breaks()) +
        facet_wrap(~"Rarefaction and accumulation") +
        theme(aspect.ratio = 1)
      
    } else {
      
      df <- data.frame("means" = c(0,mean_vec),
                       "lbs" = c(0,ci_mat[,1]),
                       "ubs" = c(0,ci_mat[,2]))
      
      p <- ggplot(df, aes(x = 0:n, y = means)) +
        geom_line(linetype = "dashed", color = "red" , linewidth = 0.9) +
        geom_ribbon(aes(ymin = lbs, ymax = ubs), color = "red", linewidth = 0.8, alpha = 0.1) +
        xlab("# observations") + ylab("# distinct features") + 
        theme_light() + 
        theme(legend.position = "top") +
        scale_y_continuous(breaks = pretty_breaks()) +
        scale_x_continuous(breaks = pretty_breaks()) +
        facet_wrap(~"Rarefaction") +
        theme(aspect.ratio = 1)
      
    }
    
    
  } else if (type == "extrapolation"){
    
    if(is.null(M)){
      M = n
    }
    
    extr_list <- extrapolation(object = x, M = M, seed = seed)
    
    ci_mat <- as.matrix(bind_rows(lapply(extr_list, quantile, probs = c(0.025,0.975) )))
    mean_vec <- as.vector(unlist(lapply(extr_list, mean)))
    
    # Extract accumulation curve of the observed sample (or average accumulation)
    accum <- rarefaction(x$feature_matrix)
    accum_df <- data.frame("accum" = c(0, accum))
    
    df <- data.frame("means" = c(Kn,mean_vec),
                     "lbs" = c(Kn,ci_mat[,1]),
                     "ubs" = c(Kn,ci_mat[,2]))
    
    p <- ggplot(df, aes(x = n:(n+M), y = means)) +
      geom_line(linetype = "dashed", color = "red" , linewidth = 0.9) +
      geom_ribbon(aes(ymin = lbs, ymax = ubs), color = "red" , linewidth = 0.8, alpha = 0.1) +
      geom_point( data = accum_df, aes(x = 0:n, y = accum), color="black", shape = 1) +
      geom_vline(aes(xintercept = n), linetype="dashed", color = "grey") +
      xlab("# observations") + ylab("# distinct features") + 
      theme_light() + 
      theme(legend.position = "top") +
      scale_y_continuous(breaks = pretty_breaks()) +
      scale_x_continuous(breaks = pretty_breaks()) +
      facet_wrap(~"Extrapolation and accumulation") +
      theme(aspect.ratio = 1)
    
  
    
  }
  
  return(p)
  
}