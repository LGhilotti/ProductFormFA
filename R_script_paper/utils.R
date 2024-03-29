
# Functions to fit PoissonBB, NegBinBB and GammaIBP on the specific training set and Nbar --------

fit_PoissonBB_wrap <- function(feature_matrix, Nbar){
  
  n_train <- nrow(feature_matrix)
  
  # Initialization and MCMC setting 
  init_PoissonBB <- list(alpha_0 = -1, s_0 = 1)
  init_obj_PoissonBB <- initialization(model = "PoissonBB", init = init_PoissonBB )
  mcmcparams_PoissonBB <- list(tau = 0.1, S = 300, n_burnin = 100, thin = 2)
  mcmcparams_obj_PoissonBB <- mcmcparameters(model = "PoissonBB", mcmcparams = mcmcparams_PoissonBB)
  
  # Hyperparameters elicitation 
  hyper_PoissonBB <- list(a_alpha = 1, b_alpha = 0.1,
                          a_s = 2, b_s = 0.2,
                          lambda = Nbar)
  prior_obj_PoissonBB <- prior(model = "PoissonBB", hyper = hyper_PoissonBB) 
  
  # Fit the model
  PoissonBB_fit <- GibbsFA(feature_matrix = feature_matrix, 
                           model = "PoissonBB", 
                           prior = prior_obj_PoissonBB, 
                           initialization = init_obj_PoissonBB, 
                           mcmcparams = mcmcparams_obj_PoissonBB)
  
  return(PoissonBB_fit)
  
}



fit_NegBinBB_wrap <- function(feature_matrix, Nbar){
  
  n_train <- nrow(feature_matrix)
  
  # Initialization and MCMC setting
  init_NegBinBB <- list(alpha_0 = -1, s_0 = 1)
  init_obj_NegBinBB <- initialization(model = "NegBinBB", init = init_NegBinBB )
  mcmcparams_NegBinBB <- list(tau = 0.1, S = 300, n_burnin = 100, thin = 2)
  mcmcparams_obj_NegBinBB <- mcmcparameters(model = "NegBinBB", mcmcparams = mcmcparams_NegBinBB)
  
  # Hyperparameters elicitation 
  c_fr <- 10
  
  hyper_NegBinBB <- list(a_alpha = 1, b_alpha = 0.1,
                         a_s = 2, b_s = 0.2,
                         n0 = Nbar/(c_fr - 1), # n0, mu0 are set s.t. E(N) = Nbar, Var(N) = c_fr*E(N)
                         mu0 = 1/c_fr)
  prior_obj_NegBinBB <- prior(model = "NegBinBB", hyper = hyper_NegBinBB) 
  
  
  # Fit the model
  NegBinBB_fit <- GibbsFA(feature_matrix = feature_matrix, 
                          model = "NegBinBB", 
                          prior = prior_obj_NegBinBB, 
                          initialization = init_obj_NegBinBB, 
                          mcmcparams = mcmcparams_obj_NegBinBB)
  
  return(NegBinBB_fit)
  
}


fit_GammaIBP_wrap <- function(feature_matrix){
  
  n_train <- nrow(feature_matrix)
  
  # Initialization and MCMC setting 
  init_GammaIBP <- list(alpha_0 = 0.5, s_0 = 15, a_0 = 5, b_0 = 1)
  init_obj_GammaIBP <- initialization(model = "GammaIBP", init = init_GammaIBP )
  mcmcparams_GammaIBP <- list(sigq_alpha = 0.1, sigq_s = 0.1, 
                              S = 300, n_burnin = 100, thin = 2)
  mcmcparams_obj_GammaIBP <- mcmcparameters(model = "GammaIBP", mcmcparams = mcmcparams_GammaIBP)
  
  # Hyperparameters elicitation 
  hyper_GammaIBP <- list(a_alpha = 2, b_alpha = 2,
                         a_s = 2, b_s = 0.2,
                         q = 0.05, r = 1, t = 0.1)
  prior_obj_GammaIBP <- prior(model = "GammaIBP", hyper = hyper_GammaIBP) 
  
  
  # Fit the model
  GammaIBP_fit <- GibbsFA(feature_matrix = feature_matrix, 
                          model = "GammaIBP", 
                          prior = prior_obj_GammaIBP, 
                          initialization = init_obj_GammaIBP, 
                          mcmcparams = mcmcparams_obj_GammaIBP)
  
  return(GammaIBP_fit)
  
}



# Function for accuracy metrics -----

compute_accuracy <- function(obs_n, est_n, obs_t) {
  # res <- 1/(1 + abs(obs_n - est_n))
  res <- abs(obs_n - est_n)/obs_t
  return ( res )
}

# Function for extrapolation plots -----

list_extr_GibbsFA_to_long <- function(list_extr, model){
  
  if (!model %in% c("Poisson", "NegBin", "Gamma")){
    stop("Error: invalid model type.")
  }
  
  df_extr <- tibble()
  
  for (j in 1:length(list_extr)){
    
    df_ci <- as_tibble(t(bind_rows(as.data.frame(lapply(list_extr[[j]], quantile, prob = c(0.025, 0.975))),
                                   as.data.frame(lapply(list_extr[[j]], mean))))) 
    colnames(df_ci) <- c("lbs", "ubs", "means")
    df_ci <- df_ci %>%
      add_column("t" = 1:nrow(df_ci)) 
    
    if (model == "Gamma"){
      df_ci <- df_ci %>% 
        add_column(Setting = names(list_extr[j])) %>%
        extract(Setting, c("n_train"), "n_train\\.([[:digit:]]+)") %>%
        add_column(Nbar = "Not applicable")
      
    } else {
      df_ci <- df_ci %>% 
        add_column(Setting = names(list_extr[j])) %>%
        extract(Setting, c("n_train","Nbar"), "n_train\\.([[:digit:]]+)\\:Nbar\\.([[:alnum:]]+)")
    }
    
    df_extr <- bind_rows(df_extr, df_ci)
  }
  
  df_extr$n_train <- as.integer(df_extr$n_train)
  df_extr <- df_extr %>%
    mutate( t = t + n_train ) %>%
    add_column("model" = model)
  
  return(df_extr)
  
}

list_extr_competitor_to_long <- function(list_extr, model){
  
  if (!model %in% c("GT", "Chao")){
    stop("Error: invalid model type.")
  }
  
  df_extr <- tibble()
  
  for (j in 1:length(list_extr)){
    
    if (model == "GT"){
      df_ci <- as_tibble(list_extr[[j]]) 
      df_ci <- df_ci %>%
        add_column("t" = 0:(nrow(df_ci)-1)) %>%
        add_column(Setting = names(list_extr[j])) %>%
        extract(Setting, c("n_train"), "n_train\\.([[:digit:]]+)") %>%
        add_column(Nbar = "Not applicable") 
    }
    
    if (model == "Chao"){
      df_ci <- as_tibble(list_extr[[j]]) 
      df_ci <- df_ci %>%
        select(-c(lbs,ubs)) %>%
        add_column(Setting = names(list_extr[j])) %>%
        extract(Setting, c("n_train"), "n_train\\.([[:digit:]]+)") %>%
        add_column(Nbar = "Not applicable") 
      
    } 
    
    df_extr <- bind_rows(df_extr, df_ci)
  }
  
  df_extr$n_train <- as.integer(df_extr$n_train)
  df_extr <- df_extr %>%
    add_column("model" = model) %>%
    filter(t >= n_train)        
  
  return(df_extr)
  
}



# Function for Beta-Binomial estimator ----

beta_binomial_estimator <- function(data_mat){
  
  # Compute total number of sites
  n <- nrow(data_mat)
  
  # Delete NA
  data_mat <- data_mat[, colSums(is.na(data_mat))==0]
  
  # Delete zero-columns
  data_mat <- data_mat[, colSums(data_mat)!=0]
  
  # Set K to be the observed number of features
  K <- ncol(data_mat)
  
  # Compute number of species contained in exactly k individuals
  counts <- colSums(data_mat)
  
  Q_1 <- sum(counts == 1)
  Q_2 <- sum(counts == 2)
  Q_3 <- sum(counts == 3)
  
  # Compute Q_hat_0 
  if (Q_2 == 0){
    Q_hat_0 <- (n-1)/n * (Q_1*(Q_1 -1))/2
  } else {
    Q_hat_0 <- (n-1)/n * (Q_1^2)/(2*Q_2)
  }
  
  if (Q_1 == 0) { Q_1 <- 1}
  if (Q_3 == 0) { Q_3 <- 1}
  
  # Compute last term
  if (2*Q_2^2 / (3*Q_1*Q_3) <= 1){
    last <- 2 - max(0.5, 2*Q_2^2 / (3*Q_1*Q_3))
  } else {
    last <- 1
  }
  
  # Compute the statistic
  res <- K + Q_hat_0 * last
  
  return (res)
  
}

# Functions for smoothed Good-Toulmin extrapolation ----

predict_good_toulmin <- function(N, M, sfs, cts, alternative = 0){
  
  preds <- rep(0,N+M+1)
  vars_ <- rep(0,N+M+1)
  preds[1:(N+1)] <- cts[1:(N+1)]
  preds_vars <- lapply(1:M, function(m) missed_gt(N, m, sfs, alternative))
  
  preds[(N+2):length(preds)] = cts[N+1] + sapply(preds_vars, function(p) p[1])
  vars_[(N+2):length(vars_)] = sapply(preds_vars, function(p) p[2])
  
  return (list("preds" = preds, "vars" = vars_))
  
}


missed_gt <- function(N, M, sfs, alternative = 0){
  
  if (length(sfs)>N){
    stop('Too many entries in the sfs; 1-th entry should be # things observed once; last entry # things observed N times')
  }
  
  signed_sfs = (-1)^(2:(length(sfs)+1)) * sfs
  t = M/N
  t_power = t^(1:length(sfs))
  if (M <= N){
    preds = sum(signed_sfs*t_power)
    vars_ = sum(sfs*(t_power^2))
  } else {
    if (alternative == T){
      kappa = floor(0.5 * log(N * (t^2) /(t-1), base = 2))
      theta = 1/(t+1)
    } else {
      kappa = floor(0.5 * log(N * (t^2) /(t-1), base = 2)/log(3))
      theta = 2/(t+1)
    }
    prob = 1-pbinom(size=kappa, prob=theta, q=0:(length(sfs)-1))
    preds = sum(signed_sfs*t_power*prob)
    vars_ = sum(abs(signed_sfs)*(t_power^2)*(prob^2))
  }
  
  
  return (c(preds, vars_))
  
}



