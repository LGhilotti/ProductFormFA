#' @export
#'
mcmcparameters <- function(model, mcmcparams){
  
  if (model == "PoissonBB" | model == "NegBinBB") {
    
    if (!all(names(mcmcparams) == c("tau","S", "n_burnin", "thin")) ){
      stop(paste0("Incorrect set of MCMC parameters for ", model))
    }
    
    if (mcmcparams$tau <= 0 | mcmcparams$n_burnin < 0 | mcmcparams$S <= mcmcparams$n_burnin | mcmcparams$thin <= 0){
      stop(paste0("Invalid value of some MCMC parameters for ", model))
    }
    
    class(mcmcparams) <- c("mcmcparameters", model)
    return(mcmcparams) 
  }
  
  
  if (model == "GammaIBP") {
    
    if (!all(names(mcmcparams) == c("sigq_alpha","sigq_s","S", "n_burnin", "thin")) ){
      stop("Incorrect set of MCMC parameters for GammaIBP.")
    }
    
    if (mcmcparams$sigq_alpha <= 0 | mcmcparams$sigq_s <= 0 | mcmcparams$n_burnin < 0 |
        mcmcparams$S <= mcmcparams$n_burnin | mcmcparams$thin <= 0){
      stop("Invalid value of some MCMC parameters for GammaIBP.")
    }
    
    class(mcmcparams) <- c("mcmcparameters", "GammaIBP")
    return(mcmcparams)
    
  }
  
}


