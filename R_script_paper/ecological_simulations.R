
# Function to generate binary matrix according to 4 different mechanisms

generate_data <- function(mechanism, n, H, seed = 1234){

  if (!mechanism %in% c("homogeneous", "random uniform", "broken stick", "log-normal")){
    stop("Invalid generating mechanism.")
  }
  
  set.seed(seed)
  
  if (mechanism == "homogeneous"){ 
    
    pres_prob <- 0.05
    
    data_mat <- matrix(rbinom(n*H, size = 1, prob = pres_prob), nrow = n, ncol = H)
    
  } else if (mechanism == "random uniform"){
    
    as <- runif(H, 0, 1)
    c <- 0.5 / max(as)
    pis <- c*as
    
    data_mat <- matrix(rbinom(n*H, size = 1, prob = rep(pis, n)),
                       nrow = n, ncol = H, byrow = T )
    
  } else if (mechanism == "broken stick"){
    
    as <- rexp(H, 1)
    c <- 0.5 / max(as)
    pis <- c*as
    
    data_mat <- matrix(rbinom(n*H, size = 1, prob = rep(pis, n)),
                       nrow = n, ncol = H, byrow = T )
    
  } else if (mechanism == "log-normal"){
    
    as <- exp(rnorm(H, 0, 1))
    c <- 1 / max(as)
    pis <- c*as
    
    data_mat <- matrix(rbinom(n*H, size = 1, prob = rep(pis, n)),
                       nrow = n, ncol = H, byrow = T )
    
  }
  
  return(data_mat)
  
}


# Function to produce simulation on specific ecological scenario, for single dataset

ecological_scenario_singledataset <- function(mechanism, seed = 1234){
  
  if (!mechanism %in% c("homogeneous", "random uniform", "broken stick", "log-normal")){
    stop("Invalid generating mechanism.")
  }
  
  # set number of individuals
  Ns <- c(20, 40, 80)
  L <- 150
  
  # Set desired value of E[N] = Nbar
  Nbars <- c(200, 400, 600)
  # Set value of c_fr such that Var(N)= c_fr * Nbar, when two or more parameters
  c_fr <- 10
  
  # Generate data -------
  data_mat <- generate_data(mechanism = mechanism, seed = seed)
  
  # PoissonBB model: initialization and MCMC setting -----
  init_PoissonBB <- list(alpha_0 = -1, s_0 = 1)
  init_obj_PoissonBB <- initialization(model = "PoissonBB", init = init_PoissonBB )
  mcmcparams_PoissonBB <- list(tau = 0.1, S = 200, n_burnin = 10, thin = 2)
  mcmcparams_obj_PoissonBB <- mcmcparameters(model = "PoissonBB", mcmcparams = mcmcparams_PoissonBB)
  
  # NegBinBB model: initialization and MCMC setting -----
  init_NegBinBB <- list(alpha_0 = -1, s_0 = 1)
  init_obj_NegBinBB <- initialization(model = "NegBinBB", init = init_NegBinBB )
  mcmcparams_NegBinBB <- list(tau = 0.1, S = 100, n_burnin = 10, thin = 2)
  mcmcparams_obj_NegBinBB <- mcmcparameters(model = "NegBinBB", mcmcparams = mcmcparams_NegBinBB)
  
  # GammaIBP model: initialization and MCMC setting -----
  init_GammaIBP <- list(alpha_0 = 0.5, s_0 = 15, a_0 = 5, b_0 = 1)
  init_obj_GammaIBP <- initialization(model = "GammaIBP", init = init_GammaIBP )
  mcmcparams_GammaIBP <- list(sigq_alpha = 0.1, sigq_s = 0.1, 
                     S = 100, n_burnin = 10, thin = 2)
  mcmcparams_obj_GammaIBP <- mcmcparameters(model = "GammaIBP", mcmcparams = mcmcparams_GammaIBP)
  
  
  # at the end, save the workspace related to the scenario just performed
  
}