#' @export
#'
initialization <- function(model, init){
  
  if (model == "PoissonBB") {
    
    if (!all(names(init) == c("alpha_0","s_0")) ){
      stop("Incorrect set of initial parameters for PoissonBB.")
    }
    
    if (init$alpha_0 >= 0 | init$s_0 <=0 ){
      stop("Invalid value of some initial parameters for PoissonBB.")
    }
    
    init["alpha_bar_0"] = - init$alpha_0
    
    class(init) <- c("initialization", "PoissonBB")
    return(init) 
  }
  
  if (model == "NegBinBB") {
    
    if (!all(names(init) == c("alpha_0","s_0")) ){
      stop("Incorrect set of initial parameters for NegBinBB.")
    }
    
    if (init$alpha_0 >= 0 | init$s_0 <=0){
      stop("Invalid value of some initial parameters for NegBinBB.")
    }
    
    init["alpha_bar_0"] = - init$alpha_0
    
    class(init) <- c("initialization", "NegBinBB")
    return(init)
    
  }
  
  if (model == "GammaIBP") {
    
    if (!all(names(init) == c("alpha_0","s_0","a_0","b_0")) ){
      stop("Incorrect set of initial parameters for GammaIBP.")
    }
    
    if (init$alpha_0 <= 0 | init$alpha_0 >= 1 | init$s_0 <=0 | init$a_0 <= 0 | init$b_0 <= 0 ){
      stop("Invalid value of some initial parameters for GammaIBP.")
    }
    
    class(init) <- c("initialization", "GammaIBP")
    return(init)
    
  }
  
}


