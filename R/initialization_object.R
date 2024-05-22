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


#' @import purrr
#' @export
#'
eb_params <- function(model, init, known){
  
  if (length(intersect(names(init), names(known))) != 0){
    stop("Some parameters are both in init and known list.")
  }
  
  if (model == "BB") {
    
    if (!setequal(union(names(init), names(known)), c("alpha","s", "Nhat_prime")) ){
      stop("Some parameters are missing for BB.")
    }
    
    null_known <- rep(NA, length(names(init)))
    names(null_known) <- names(init)
    null_known <- c(as_vector(known), null_known)
    
    init <- c(as_vector(init), as_vector(known) ) # here I have both init and known values
    
    if (init["alpha"] >= 0 | init["s"] <=0 | init["Nhat_prime"] <= 0){
      stop("Invalid value of some initial parameters for BB.")
    }
    
    # order the list as (alpha, s, Nhat_prime)
    par <- list("init" = init[order(factor(names(init), 
                                           levels=c("alpha", "s", "Nhat_prime")))],
                "known" = null_known[order(factor(names(null_known), 
                                                  levels=c("alpha", "s", "Nhat_prime")))])
    
    class(par) <- c("eb_params", "BB")
    return(par) 
  }
  
  if (model == "IBP") {
    
    if (!setequal(union(names(init), names(known)), c("alpha","s", "Gamma")) ){
      stop("Some parameters are missing for IBP.")
    }
    
    null_known <- rep(NA, length(names(init)))
    names(null_known) <- names(init)
    null_known <- c(as_vector(known), null_known)
    
    init <- c(as_vector(init), as_vector(known) ) # here I have both init and known values
    
    if (init["alpha"] < 0 | init["alpha"] > 1 | init["s"] <=0 | init["Gamma"] <= 0){
      stop("Invalid value of some initial parameters for IBP.")
    }
    
    # order the list as (alpha, s, Gamma)
    par <- list("init" = init[order(factor(names(init), 
                                           levels=c("alpha", "s", "Gamma")))],
                "known" = null_known[order(factor(names(null_known), 
                                                  levels=c("alpha", "s", "Gamma")))])
    
    class(par) <- c("eb_params", "IBP")
    return(par) 
  }
  
  
  if (model == "PoissonBB") {
    
    if (!setequal(union(names(init), names(known)), c("alpha","s", "lambda")) ){
      stop("Some parameters are missing for PoissonBB.")
    }
  
    null_known <- rep(NA, length(names(init)))
    names(null_known) <- names(init)
    null_known <- c(as_vector(known), null_known)
    
    init <- c(as_vector(init), as_vector(known) ) # here I have both init and known values
    
    if (init["alpha"] >= 0 | init["s"] <=0 | init["lambda"] <= 0){
      stop("Invalid value of some initial parameters for PoissonBB.")
    }
    
    # order the list as (alpha, s, lambda)
    par <- list("init" = init[order(factor(names(init), 
                                           levels=c("alpha", "s", "lambda")))],
                "known" = null_known[order(factor(names(null_known), 
                                                  levels=c("alpha", "s", "lambda")))])
    
    class(par) <- c("eb_params", "PoissonBB")
    return(par) 
  }
  
  if (model == "NegBinBB") {
    
    if (!setequal(union(names(init), names(known)), c("alpha","s", "var_fct", "mu0")) ){
      stop("Some parameters are missing for NegBinBB.")
    }
    
    null_known <- rep(NA, length(names(init)))
    names(null_known) <- names(init)
    null_known <- c(as_vector(known), null_known)
    
    init <- c(as_vector(init), as_vector(known) ) # here I have both init and known values
    
    if (init["alpha"] >= 0 | init["s"] <=0 | init["var_fct"] <= 1 | 
        init["mu0"] <= 0 ){
      stop("Invalid value of some initial parameters for NegBinBB.")
    }
    
    # order the list as (alpha, s, n0, mu0)
    par <- list("init" = init[order(factor(names(init), 
                                           levels=c("alpha", "s", "var_fct", "mu0")))],
                "known" = null_known[order(factor(names(null_known), 
                                                  levels=c("alpha", "s", "var_fct", "mu0")))])
    
    class(par) <- c("eb_params", "NegBinBB")
    return(par) 
    
  }
  
  if (model == "NegBinBB_np") {
    
    if (!setequal(union(names(init), names(known)), c("alpha","s", "n0", "p")) ){
      stop("Some parameters are missing for NegBinBB.")
    }
    
    null_known <- rep(NA, length(names(init)))
    names(null_known) <- names(init)
    null_known <- c(as_vector(known), null_known)
    
    init <- c(as_vector(init), as_vector(known) ) # here I have both init and known values
    
    if (init["alpha"] >= 0 | init["s"] <=0 | init["n0"] <= 0 | 
        init["p"] <= 0 | init["p"] >= 1){
      stop("Invalid value of some initial parameters for NegBinBB.")
    }
    
    # order the list as (alpha, s, n0, p)
    par <- list("init" = init[order(factor(names(init), 
                                           levels=c("alpha", "s", "n0", "p")))],
                "known" = null_known[order(factor(names(null_known), 
                                                  levels=c("alpha", "s", "n0", "p")))])
    
    class(par) <- c("eb_params", "NegBinBB_np")
    return(par) 
    
  }
  
  if (model == "GammaIBP") {
    
    if (!setequal(union(names(init), names(known)), c("alpha","s", "a", "b")) ){
      stop("Some parameters are missing for GammaIBP.")
    }
    
    null_known <- rep(NA, length(names(init)))
    names(null_known) <- names(init)
    null_known <- c(as_vector(known), null_known)
    
    init <- c(as_vector(init), as_vector(known) ) # here I have both init and known values
    
    if (init["alpha"] <= 0 | init["alpha"] >= 1 | init["s"] <=0 | 
        init["a"] <= 0 | init["b"] <= 0 ){
      stop("Invalid value of some initial parameters for GammaIBP.")
    }
    
    # order the list as (alpha, s, a, b)
    par <- list("init" = init[order(factor(names(init), 
                                           levels=c("alpha", "s", "a", "b")))],
                "known" = null_known[order(factor(names(null_known), 
                                                  levels=c("alpha", "s", "a", "b")))])
    
    class(par) <- c("eb_params", "GammaIBP")
    return(par) 
    
  }
  
}


