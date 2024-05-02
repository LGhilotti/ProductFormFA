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
  
  if (model == "GammaIBP") {
    
    if (!all(names(hyper) == c("a_alpha","b_alpha","a_s","b_s","q","r","t")) ){
      stop("Incorrect set of hyperparameters for GammaIBP.")
    }
    
    if (hyper$a_alpha <= 0 | hyper$b_alpha <=0 | hyper$a_s <= 0 | hyper$b_s <= 0 | 
        hyper$q <= 0 | hyper$q >= 1 | hyper$r <=0 | hyper$t <= 0){
      stop("Invalid value of some hyperparameters for GammaIBP.")
    }
    
    class(hyper) <- c("prior", "GammaIBP")
    return(hyper)
    
  }
  
}


#' @export
#'
summary.PoissonBB <- function(object, ...) {
  
  
  cat("Model:",
      "\n\t BB with Poisson(lambda) mixture",
      "\nPrior details:\n",
      paste0("\t E(alpha) = ", - object$a_alpha/ object$b_alpha, "; "),
      paste0("Var(alpha) = ", object$a_alpha/ (object$b_alpha^2)),
      paste0("\n\t E(s) = ", object$a_s/ object$b_s, "; "),
      paste0("Var(s) = ", object$a_s/ (object$b_s^2)),
      paste0("\n\t E(N) = ", object$lambda, "; "),
      paste0("Var(N) = ", object$lambda)
  )
  
}


#' @export
#'
summary.NegBinBB <- function(object, ...) {
  
  
  cat("Model:",
      "\n\t BB with NB(n0,mu0) mixture",
      "\nPrior details:\n",
      paste0("\t E(alpha) = ", - object$a_alpha/ object$b_alpha, "; "),
      paste0("Var(alpha) = ", object$a_alpha/ (object$b_alpha^2)),
      paste0("\n\t E(s) = ", object$a_s/ object$b_s, "; "),
      paste0("Var(s) = ", object$a_s/ (object$b_s^2)),
      paste0("\n\t E(N) = ", object$mu0, "; "),
      paste0("Var(N) = ", object$mu0 + object$mu0^2/ object$n0)
  )
  
}


#' @export
#'
summary.GammaIBP <- function(object, ...) {
  
  
  cat("Model:",
      "\n\t IBP with Gamma(a,b) mixture",
      "\nPrior details:\n",
      paste0("\t E(alpha) = ", object$a_alpha/ (object$a_alpha + object$b_alpha), "; "),
      paste0("Var(alpha) = ", (object$a_alpha*object$b_alpha)/ (object$a_alpha + object$b_alpha)^2 /
               (object$a_alpha + object$b_alpha + 1)),
      paste0("\n\t E(s) = ", object$a_s/ object$b_s, "; "),
      paste0("Var(s) = ", object$a_s/ (object$b_s^2)),
      paste0("\n\t E(a) = ", 1/ object$q , "; "),
      paste0("Var(a) = ", (1-object$q)/ (object$q^2) ),
      paste0("\n\t E(b) = ", object$r/ object$t , "; "),
      paste0("Var(b) = ", object$r/ object$t^2 )
  )
  
}