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
            dbeta(alpha, data_full$a, data_full$b, log = TRUE) +
            dgamma(s, data_full$a_s, data_full$b_s, log = TRUE) )
  
}


#' Compute Bayes Factor (BF) for: NegBinBB vs GammaIBP (Fully-Bayesian approach)
#'
#' @param NegBinBB_fit object of class \code{GibbsFA, NegBinBB}
#' @param GammaIBP_fit object of class \code{GibbsFA, GammaIBP}
#'
#' @export
#'
compute_BF <- function(NegBinBB_fit, GammaIBP_fit){
  
  ### Let assume that NegBinBB = H0 and GammaIBP = H1
  
  data_summary <- vector("list")
  # 1) Check if models refer to same data
  if (all(dim(NegBinBB_fit$feature_matrix) == dim(GammaIBP_fit$feature_matrix)) &
      all(NegBinBB_fit$feature_matrix == GammaIBP_fit$feature_matrix)){
    
    Z <- GammaIBP_fit$feature_matrix[, colSums(is.na(GammaIBP_fit$feature_matrix))==0]
    Z <- Z[, colSums(Z)!=0]
    counts <- colSums(Z)
    data_summary[["counts"]] <- counts
    data_summary[["n"]] <- nrow(Z)
    data_summary[["K"]] <- ncol(Z)
    
  } else{
    stop("Models are trained on different datasets")
  }
  
  # 2) Posterior samples and data with hyperparameters
  
  # NegBinBB: use parametrization (alpha_bar, s)
  samples_NegBinBB_list <- list("alpha_bar" = - NegBinBB_fit$alpha_chain, 
                           "s" = NegBinBB_fit$alpha_chain + NegBinBB_fit$theta_chain)
  
  samples_NegBinBB <- as.matrix(as.data.frame(samples_NegBinBB_list))
  
  data_full_NegBinBB <- append(data_summary,
        list("n0" = NegBinBB_fit$prior$n0,
         "mu0" = NegBinBB_fit$prior$mu0,
         "a_alpha" = NegBinBB_fit$prior$a_alpha,
         "b_alpha" = NegBinBB_fit$prior$b_alpha,
         "a_s" = NegBinBB_fit$prior$a_s,
         "b_s" = NegBinBB_fit$prior$b_s))
  
  # GammaIBP: use parametrization (alpha, s)
  samples_GammaIBP_list <- list("alpha" =  GammaIBP_fit$alpha_chain, 
                                "s" = GammaIBP_fit$alpha_chain + GammaIBP_fit$theta_chain)
  
  samples_GammaIBP <- as.matrix(as.data.frame(samples_GammaIBP_list))  
  
  data_full_GammaIBP <- append(data_summary,
                               list("a" = GammaIBP_fit$prior$a,
                                    "b" = GammaIBP_fit$prior$b,
                                    "a_alpha" = GammaIBP_fit$prior$a_alpha,
                                    "b_alpha" = GammaIBP_fit$prior$b_alpha,
                                    "a_s" = GammaIBP_fit$prior$a_s,
                                    "b_s" = GammaIBP_fit$prior$b_s))
  
  # 3) Specify parameter bounds 
  # NegBinBB
  cn <- colnames(samples_NegBinBB)
  lb_NegBinBB <- c(0, 0)
  ub_NegBinBB <- c(Inf, Inf)
  names(lb_NegBinBB) <- names(ub_NegBinBB) <- cn

  # GammaIBP
  cn <- colnames(samples_GammaIBP)
  lb_GammaIBP <- c(0, 0)
  ub_GammaIBP <- c(1, Inf)
  names(lb_GammaIBP) <- names(ub_GammaIBP) <- cn
  
  
  # 4) Compute log marginal likelihood via bridge sampling 
  # NegBinBB
  NegBinBB.bridge <- bridge_sampler(samples = samples_NegBinBB, data = data_full_NegBinBB,
                              log_posterior = log_posterior_NegBinBB_cpp, 
                              lb = lb_NegBinBB,
                              ub = ub_NegBinBB, silent = TRUE)
  print(NegBinBB.bridge)
  
  # GammaIBP
  GammaIBP.bridge <- bridge_sampler(samples = samples_GammaIBP, data = data_full_GammaIBP,
                              log_posterior = log_posterior_GammaIBP_cpp, lb = lb_GammaIBP,
                              ub = ub_GammaIBP, silent = TRUE)
  print(GammaIBP.bridge)
  
  # compute percentage error
  print(error_measures(NegBinBB.bridge)$percentage)
  print(error_measures(GammaIBP.bridge)$percentage)
  
  # compute Bayes factor
  BF01 <- bf(NegBinBB.bridge, GammaIBP.bridge, log = T)
  print(BF01)
  
  return(BF01)
  
}



#' Compute AIC/BIC: available models PoissonBB, NegBinBB, GammaIBP (Empirical-Bayes approach)
#'
#' @param eb_model_fit object of class \code{GibbsFA, PoissonBB_eb} 
#' or \code{GibbsFA, NegBinBB_eb} or \code{GibbsFA, GammaIBP_eb}
#'
#' @export
#'
compute_AICs_BICs <- function(eb_model_fit){
  
  if (!all(class(eb_model_fit) == c("GibbsFA", "PoissonBB_eb")) &
      !all(class(eb_model_fit) == c("GibbsFA", "NegBinBB_eb")) &
      !all(class(eb_model_fit) == c("GibbsFA", "GammaIBP_eb"))){
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
  
  return(list("AIC" = AIC,
              "BIC" = BIC))
  
}



