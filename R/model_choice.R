#
# Functions for model choice between NegBinBB vs GammaIBP
#


#' Log-evaluation of un-normalized NegBinBB posterior
#'
#' @param post_sample contains posterior samples
#' @param data_full contains data and all hyperparameters
#'
#'
log_posterior_NegBinBB_cpp <- function(post_sample, data_full) {
  
  alpha_bar <- post_sample["alpha_bar"]
  s <- post_sample["s"]
  
  pars <- c("alpha"= - alpha_bar, "s" = s, 
           "n0" = data_full$n0, "p" = data_full$n0 / (data_full$n0 + data_full$mu0))
  
  return ( - neg_log_EFPF_NegBinBB(n = data_full$n, counts = data_full$counts,
                            pars = pars) +
            dgamma(alpha_bar, shape = data_full$a_alpha, rate = data_full$b_alpha, log = TRUE) +
            dgamma(s, shape = data_full$a_s, rate = data_full$b_s, log = TRUE) )
  
}

#' Log-evaluation of un-normalized classicBB posterior
#'
#' @param post_sample contains posterior samples
#' @param data_full contains data and all hyperparameters
#'
#'
log_posterior_classicBB_cpp <- function(post_sample, data_full) {
  
  alpha_bar <- post_sample["alpha_bar"]
  s <- post_sample["s"]
  
  pars <- c("alpha"= - alpha_bar, "s" = s, 
            "Nhat_prime" = data_full$N - data_full$data_summary$K)
  
  return ( - neg_log_EFPF_BB(n = data_full$n, counts = data_full$counts,
                                   pars = pars) +
             dgamma(alpha_bar, shape = data_full$a_alpha, rate = data_full$b_alpha, log = TRUE) +
             dgamma(s, shape = data_full$a_s, rate = data_full$b_s, log = TRUE) )
  
}

#' Log-evaluation of un-normalized PoissonBB posterior
#'
#' @param post_sample contains posterior samples
#' @param data_full contains data and all hyperparameters
#'
#'
log_posterior_PoissonBB_cpp <- function(post_sample, data_full) {
  
  alpha_bar <- post_sample["alpha_bar"]
  s <- post_sample["s"]
  
  pars <- c("alpha"= - alpha_bar, "s" = s, 
            "lambda" = data_full$lambda)
  
  return ( - neg_log_EFPF_PoissonBB(n = data_full$n, counts = data_full$counts,
                             pars = pars) +
             dgamma(alpha_bar, shape = data_full$a_alpha, rate = data_full$b_alpha, log = TRUE) +
             dgamma(s, shape = data_full$a_s, rate = data_full$b_s, log = TRUE) )
  
}

#' Log-evaluation of un-normalized GammaIBP posterior
#'
#' @param post_sample contains posterior samples
#' @param data_full contains data and all hyperparameters
#'
#'
log_posterior_GammaIBP_cpp <- function(post_sample, data_full) {
  
  alpha <- post_sample["alpha"]
  s <- post_sample["s"]
  
  pars <- c("alpha"= alpha, "s" = s, 
            "a" = data_full$a, "b" = data_full$b)
  
  return ( - neg_log_EFPF_GammaIBP(n = data_full$n, counts = data_full$counts,
                                   pars = pars) +
            dbeta(alpha, data_full$a_alpha, data_full$b_alpha, log = TRUE) +
            dgamma(s, shape = data_full$a_s, rate = data_full$b_s, log = TRUE) )
  
}


#' Log-evaluation of un-normalized classicIBP posterior
#'
#' @param post_sample contains posterior samples
#' @param data_full contains data and all hyperparameters
#'
#'
log_posterior_classicIBP_cpp <- function(post_sample, data_full) {
  
  alpha <- post_sample["alpha"]
  s <- post_sample["s"]
  
  pars <- c("alpha"= alpha, "s" = s, 
            "Gam" = data_full$gam)
  
  return ( - neg_log_EFPF_IBP(n = data_full$n, counts = data_full$counts,
                                   pars = pars) +
             dbeta(alpha, data_full$a_alpha, data_full$b_alpha, log = TRUE) +
             dgamma(s, shape = data_full$a_s, rate = data_full$b_s, log = TRUE) )
  
}



#' Computes log marginal likelihood via bridge sampling (Fully-Bayesian approach)
#'
#' @param model_fit object of class \code{classicBB} (BB with N total features),
#' \code{PoissonBB} (BB with Poisson(lambda) mixture),
#' \code{NegBinBB} (BB with NB(n0, mu0) mixture), 
#' \code{classicIBP} (IBP with gam total mass) and \code{GammaIBP} (IBP with Gamma(a, b) mixture)
#'
#' @export
#' @import bridgesampling
#'
compute_log_marginal_likelihood_bridge <- function(model_fit){
  
  
  if (!all(class(model_fit) == c("GibbsFA", "NegBinBB")) &
      !all(class(model_fit) == c("GibbsFA", "classicBB")) &
      !all(class(model_fit) == c("GibbsFA", "PoissonBB")) &
      !all(class(model_fit) == c("GibbsFA", "GammaIBP")) &
      !all(class(model_fit) == c("GibbsFA", "classicIBP"))){
    stop("Incompatible class object")
  }
  
  data_summary <- vector("list")
  Z <- model_fit$feature_matrix[, colSums(is.na(model_fit$feature_matrix))==0]
  Z <- Z[, colSums(Z)!=0]
  counts <- colSums(Z)
  data_summary[["counts"]] <- counts
  data_summary[["n"]] <- nrow(Z)
  data_summary[["K"]] <- ncol(Z)
  
  
  # NegBinBB: use parametrization (alpha_bar, s)
  if (class(model_fit)[2] == "NegBinBB"){
    
    # 2) Posterior samples and data with hyperparameters
    samples_NegBinBB_list <- list("alpha_bar" = - model_fit$alpha_chain, 
                                  "s" = model_fit$alpha_chain + model_fit$theta_chain)
    
    samples_NegBinBB <- as.matrix(as.data.frame(samples_NegBinBB_list))
    
    data_full_NegBinBB <- append(data_summary,
                                 list("n0" = model_fit$prior$n0,
                                      "mu0" = model_fit$prior$mu0,
                                      "a_alpha" = model_fit$prior$a_alpha,
                                      "b_alpha" = model_fit$prior$b_alpha,
                                      "a_s" = model_fit$prior$a_s,
                                      "b_s" = model_fit$prior$b_s))
    
    # 3) Specify parameter bounds 
    cn <- colnames(samples_NegBinBB)
    lb_NegBinBB <- c(0, 0)
    ub_NegBinBB <- c(Inf, Inf)
    names(lb_NegBinBB) <- names(ub_NegBinBB) <- cn
    
    # 4) Compute log marginal likelihood via bridge sampling 
    model.bridge <- bridge_sampler(samples = samples_NegBinBB, data = data_full_NegBinBB,
                                      log_posterior = log_posterior_NegBinBB_cpp, 
                                      lb = lb_NegBinBB,
                                      ub = ub_NegBinBB, silent = TRUE)
    
  }
  
  # classicBB: use parametrization (alpha_bar, s)
  if (class(model_fit)[2] == "classicBB"){
    
    # 2) Posterior samples and data with hyperparameters
    samples_classicBB_list <- list("alpha_bar" = - model_fit$alpha_chain, 
                                  "s" = model_fit$alpha_chain + model_fit$theta_chain)
    
    samples_classicBB <- as.matrix(as.data.frame(samples_classicBB_list))
    
    data_full_classicBB <- append(data_summary,
                                 list("N" = model_fit$prior$N,
                                      "a_alpha" = model_fit$prior$a_alpha,
                                      "b_alpha" = model_fit$prior$b_alpha,
                                      "a_s" = model_fit$prior$a_s,
                                      "b_s" = model_fit$prior$b_s))
    
    # 3) Specify parameter bounds 
    cn <- colnames(samples_classicBB)
    lb_classicBB <- c(0, 0)
    ub_classicBB <- c(Inf, Inf)
    names(lb_classicBB) <- names(ub_classicBB) <- cn
    
    # 4) Compute log marginal likelihood via bridge sampling 
    model.bridge <- bridge_sampler(samples = samples_classicBB, data = data_full_classicBB,
                                   log_posterior = log_posterior_classicBB_cpp, 
                                   lb = lb_classicBB,
                                   ub = ub_classicBB, silent = TRUE)
    
  }
  
  # PoissonBB: use parametrization (alpha_bar, s)
  if (class(model_fit)[2] == "PoissonBB"){
    
    # 2) Posterior samples and data with hyperparameters
    samples_PoissonBB_list <- list("alpha_bar" = - model_fit$alpha_chain, 
                                   "s" = model_fit$alpha_chain + model_fit$theta_chain)
    
    samples_PoissonBB <- as.matrix(as.data.frame(samples_PoissonBB_list))
    
    data_full_PoissonBB <- append(data_summary,
                                  list("lambda" = model_fit$prior$lambda,
                                       "a_alpha" = model_fit$prior$a_alpha,
                                       "b_alpha" = model_fit$prior$b_alpha,
                                       "a_s" = model_fit$prior$a_s,
                                       "b_s" = model_fit$prior$b_s))
    
    # 3) Specify parameter bounds 
    cn <- colnames(samples_PoissonBB)
    lb_PoissonBB <- c(0, 0)
    ub_PoissonBB <- c(Inf, Inf)
    names(lb_PoissonBB) <- names(ub_PoissonBB) <- cn
    
    # 4) Compute log marginal likelihood via bridge sampling 
    model.bridge <- bridge_sampler(samples = samples_PoissonBB, data = data_full_PoissonBB,
                                   log_posterior = log_posterior_PoissonBB_cpp, 
                                   lb = lb_PoissonBB,
                                   ub = ub_PoissonBB, silent = TRUE)
    
  }
  
  # GammaIBP: use parametrization (alpha, s)
  if (class(model_fit)[2] == "GammaIBP"){
    
    # 2) Posterior samples and data with hyperparameters
    samples_GammaIBP_list <- list("alpha" =  model_fit$alpha_chain, 
                                  "s" = model_fit$alpha_chain + model_fit$theta_chain)
    
    samples_GammaIBP <- as.matrix(as.data.frame(samples_GammaIBP_list))  
    
    data_full_GammaIBP <- append(data_summary,
                                 list("a" = model_fit$prior$a,
                                      "b" = model_fit$prior$b,
                                      "a_alpha" = model_fit$prior$a_alpha,
                                      "b_alpha" = model_fit$prior$b_alpha,
                                      "a_s" = model_fit$prior$a_s,
                                      "b_s" = model_fit$prior$b_s))
    
    # 3) Specify parameter bounds 
    cn <- colnames(samples_GammaIBP)
    lb_GammaIBP <- c(0, 0)
    ub_GammaIBP <- c(1, Inf)
    names(lb_GammaIBP) <- names(ub_GammaIBP) <- cn
    
    # 4) Compute log marginal likelihood via bridge sampling 
    model.bridge <- bridge_sampler(samples = samples_GammaIBP, data = data_full_GammaIBP,
                                      log_posterior = log_posterior_GammaIBP_cpp, lb = lb_GammaIBP,
                                      ub = ub_GammaIBP, silent = TRUE)
    
    
  }
  
  
  # classicIBP: use parametrization (alpha, s)
  if (class(model_fit)[2] == "classicIBP"){
    
    # 2) Posterior samples and data with hyperparameters
    samples_classicIBP_list <- list("alpha" =  model_fit$alpha_chain, 
                                  "s" = model_fit$alpha_chain + model_fit$theta_chain)
    
    samples_classicIBP <- as.matrix(as.data.frame(samples_classicIBP_list))  
    
    data_full_classicIBP <- append(data_summary,
                                 list("gam" = model_fit$prior$gam,
                                      "a_alpha" = model_fit$prior$a_alpha,
                                      "b_alpha" = model_fit$prior$b_alpha,
                                      "a_s" = model_fit$prior$a_s,
                                      "b_s" = model_fit$prior$b_s))
    
    # 3) Specify parameter bounds 
    cn <- colnames(samples_classicIBP)
    lb_classicIBP <- c(0, 0)
    ub_classicIBP <- c(1, Inf)
    names(lb_classicIBP) <- names(ub_classicIBP) <- cn
    
    # 4) Compute log marginal likelihood via bridge sampling 
    model.bridge <- bridge_sampler(samples = samples_classicIBP, data = data_full_classicIBP,
                                   log_posterior = log_posterior_classicIBP_cpp, lb = lb_classicIBP,
                                   ub = ub_classicIBP, silent = TRUE)
    
    
  }
  
  
  return(model.bridge)
  
}



#' Compute AIC/BIC (Empirical-Bayes approach)
#'
#' @param eb_model_fit object of class \code{classicBB_eb} (BB with N total features),
#' \code{PoissonBB_eb} (BB with Poisson(lambda) mixture),
#' \code{NegBinBB_eb} (BB with NB(n0, mu0) mixture), 
#' \code{classicIBP_eb} (IBP with gam total mass) and \code{GammaIBP_eb} (IBP with Gamma(a, b) mixture)
#'
#' @export
#'
compute_AICs_BICs <- function(eb_model_fit){
  
  if (!all(class(eb_model_fit) == c("GibbsFA", "NegBinBB_eb")) &
      !all(class(eb_model_fit) == c("GibbsFA", "classicBB_eb")) &
      !all(class(eb_model_fit) == c("GibbsFA", "PoissonBB_eb")) &
      !all(class(eb_model_fit) == c("GibbsFA", "GammaIBP_eb")) &
      !all(class(eb_model_fit) == c("GibbsFA", "classicIBP_eb"))){
    stop("Incompatible class object")
  }

  data_summary <- vector("list")
  Z <- eb_model_fit$feature_matrix[, colSums(is.na(eb_model_fit$feature_matrix))==0]
  Z <- Z[, colSums(Z)!=0]
  counts <- colSums(Z)
  data_summary[["counts"]] <- counts
  data_summary[["n"]] <- nrow(Z)
  data_summary[["K"]] <- ncol(Z)
    

  # 2) Compute AIC/BIC 
  # NegBinBB
  if (class(eb_model_fit)[2] == "NegBinBB_eb"){
    pars <- c("alpha"= eb_model_fit$alpha, 
              "s" = eb_model_fit$alpha + eb_model_fit$theta, 
              "n0" = eb_model_fit$n0, 
              "p" = eb_model_fit$n0 / (eb_model_fit$n0 + eb_model_fit$mu0))
    
    max_log_efpf <- - neg_log_EFPF_NegBinBB(n = data_summary$n,
                                            counts = data_summary$counts,
                                            pars = pars)
    
    AIC <- 2*length(pars - 1) - 2*max_log_efpf
    BIC <- log(data_summary$n)*length(pars - 1) - 2*max_log_efpf
    
  }
  # PoissonBB
  if (class(eb_model_fit)[2] == "PoissonBB_eb"){
    pars <- c("alpha"= eb_model_fit$alpha, 
              "s" = eb_model_fit$alpha + eb_model_fit$theta, 
              "lambda" = eb_model_fit$lambda)
    
    max_log_efpf <- - neg_log_EFPF_PoissonBB(n = data_summary$n,
                                             counts = data_summary$counts,
                                             pars = pars)
    
    AIC <- 2*length(pars) - 2*max_log_efpf
    BIC <- log(data_summary$n)*length(pars) - 2*max_log_efpf
    
  }
  # classicBB
  if (class(eb_model_fit)[2] == "classicBB_eb"){
    pars <- c("alpha"= eb_model_fit$alpha, 
              "s" = eb_model_fit$alpha + eb_model_fit$theta, 
              "Nhat_prime" = eb_model_fit$N - data_summary$K)
    
    max_log_efpf <- - neg_log_EFPF_BB(n = data_summary$n,
                                      counts = data_summary$counts,
                                      pars = pars)
    
    AIC <- 2*length(pars) - 2*max_log_efpf
    BIC <- log(data_summary$n)*length(pars) - 2*max_log_efpf
    
  }
  # GammaIBP
  if (class(eb_model_fit)[2] == "GammaIBP_eb"){
    pars <- c("alpha"= eb_model_fit$alpha, 
              "s" = eb_model_fit$alpha + eb_model_fit$theta, 
              "a" = eb_model_fit$a, 
              "b" = eb_model_fit$b)
    
    max_log_efpf <- - neg_log_EFPF_GammaIBP(n = data_summary$n,
                                            counts = data_summary$counts,
                                            pars = pars)
    
    AIC <- 2*length(pars - 1) - 2*max_log_efpf
    BIC <- log(data_summary$n)*length(pars - 1) - 2*max_log_efpf
    
  }
  # classicIBP
  if (class(eb_model_fit)[2] == "classicIBP_eb"){
    pars <- c("alpha"= eb_model_fit$alpha, 
              "s" = eb_model_fit$alpha + eb_model_fit$theta, 
              "gam" = eb_model_fit$gam)
    
    max_log_efpf <- - neg_log_EFPF_IBP(n = data_summary$n,
                                       counts = data_summary$counts,
                                       pars = pars)
    
    AIC <- 2*length(pars - 1) - 2*max_log_efpf
    BIC <- log(data_summary$n)*length(pars - 1) - 2*max_log_efpf
    
  }
  
  return(list("AIC" = AIC,
              "BIC" = BIC,
              "max_log_efpf" = max_log_efpf))
  
}



