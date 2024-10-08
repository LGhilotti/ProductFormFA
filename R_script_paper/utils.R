
# Function for extrapolation plots -----

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



