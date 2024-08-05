#' @export
#'
prior <- function(model, hyper){
  
  if (model == "PoissonBB") {
    
    if (!all(names(hyper) == c("a_alpha","b_alpha","a_s","b_s","lambda")) ){
      stop("Incorrect set of hyperparameters for PoissonBB.")
    }
    
    if (hyper$a_alpha <= 0 | hyper$b_alpha <=0 | hyper$a_s <= 0 | hyper$b_s <= 0 | hyper$lambda <= 0){
      stop("Invalid value of some hyperparameters for PoissonBB.")
    }
    
    class(hyper) <- c("prior", "PoissonBB")
    return(hyper) 
  }
  
  if (model == "NegBinBB") {
    
    if (!all(names(hyper) == c("a_alpha","b_alpha","a_s","b_s","n0","mu0")) ){
      stop("Incorrect set of hyperparameters for NegBinBB.")
    }
    
    if (hyper$a_alpha <= 0 | hyper$b_alpha <=0 | hyper$a_s <= 0 | hyper$b_s <= 0 | 
        hyper$n0 <= 0 | hyper$mu0 <=0 ){
      stop("Invalid value of some hyperparameters for NegBinBB.")
    }
    
    class(hyper) <- c("prior", "NegBinBB")
    return(hyper)
    
  }
  
  if (model == "GammaIBP_more_prior") {
    
    if (!all(names(hyper) == c("a_alpha","b_alpha","a_s","b_s","q","r","t")) ){
      stop("Incorrect set of hyperparameters for GammaIBP (more prior).")
    }
    
    if (hyper$a_alpha <= 0 | hyper$b_alpha <=0 | hyper$a_s <= 0 | hyper$b_s <= 0 | 
        hyper$q <= 0 | hyper$q >= 1 | hyper$r <=0 | hyper$t <= 0){
      stop("Invalid value of some hyperparameters for GammaIBP (more prior).")
    }
    
    class(hyper) <- c("prior", "GammaIBP_more_prior")
    return(hyper)
    
  }
  
  
  if (model == "GammaIBP_single_prior") {
    
    if (!all(names(hyper) == c("a", "b", "a_alpha","b_alpha","a_s","b_s")) ){
      stop("Incorrect set of hyperparameters for GammaIBP (single prior).")
    }
    
    if (hyper$a <= 0 | hyper$b <=0 | 
        hyper$a_alpha <= 0 | hyper$b_alpha <=0 | hyper$a_s <= 0 | hyper$b_s <= 0){
      stop("Invalid value of some hyperparameters for GammaIBP (single prior).")
    }
    
    class(hyper) <- c("prior", "GammaIBP_single_prior")
    return(hyper)
    
  }
  
}


