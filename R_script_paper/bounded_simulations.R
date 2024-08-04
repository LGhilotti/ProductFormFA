#
# Bounded scenario simulations ####
#

rm(list=ls())
library(ProductFormFA)
library(tidyverse)
library(scales)
library(ggpubr)
library(ggthemes)
library(patchwork)
library(latex2exp)
library(dplyr, warn.conflicts = FALSE)

source("R_script_paper/Routine_Chao.R")
source("R_script_paper/utils.R")


# Function to generate binary matrix for scenario A and C ----

generate_data <- function(mechanism, n, H = 100, seed = 1234){
  
  if (!mechanism %in% c("custom", "beta_pis")){
    stop("Invalid generating mechanism.")
  }
  
  set.seed(seed)
  
  if (mechanism == "custom"){
    
    mul <- 1
    data_mat <- cbind(
    
      matrix(rbinom(n * 200*mul, size = 1, prob = 0.015), n, 200*mul),
      matrix(rbinom(n * 200*mul, size = 1, prob = 0.01), n, 200*mul),
      matrix(rbinom(n * 200*mul, size = 1, prob = 0.005), n, 200*mul)
      
    )

  } else if (mechanism == "beta_pis"){
    
    pis <- rbeta(H, 1, 100) 
    
    data_mat <- matrix(rbinom(n*H, size = 1, prob = rep(pis, n)),
                       nrow = n, ncol = H, byrow = T )
    
  } 
  
  return(data_mat)
  
}



# Function to produce fit and estimate on specific bounded scenarios ------

eb_EFPF_fit_estimate_bounded_scenario <- function(mechanism,
                                                  vars_fct_NegBinBB,
                                                  vars_GammaIBP,
                                                  eb_init_BB, eb_known_BB,
                                                  eb_init_IBP, eb_known_IBP,
                                                  seed = 1234){
  
  if (!mechanism %in% c("custom", "beta_pis")){
    stop("Invalid generating mechanism.")
  }
  
  # Set training dimensions and dimension of the whole sample 
  if (mechanism == "custom"){
    Ns <- c(50, 250, 1250)
  } else {
    Ns <- c(200, 1000, 5000) 
  }
  N_max <- Ns[length(Ns)]
  L <- N_max + 600
  
  # Set maximum number of features
  H <- 500

  # Generate data 
  data_mat <- generate_data(mechanism = mechanism, n = L, H = H, seed = seed)
  
  if (mechanism == "custom"){
    H <- ncol(data_mat) # total features, counting absent ones
  }
  
  data_mat <- data_mat[, colSums(data_mat) > 0]
  
  Nbars <- c(200, 400, 800)
  # Set grid of desired value of E[N] = Nbar, with the empirical case to be added
  if (mechanism == "custom"){
    H_hundred <- round(H / 100) * 100
    Nbars <- c(H_hundred - 400, H_hundred - 200, H_hundred + 200)
  }
  
  # Structures to save fit and estimates 
  
  # List to store the fitted models for PoissonBB, NegBinBB and GammaIBP, under different settings 
  list_eb_EFPF_fit_PoissonBB <- vector(mode = "list")
  
  list_eb_EFPF_fit_NegBinBB <-  vector(mode = "list", length = length(vars_fct_NegBinBB))
  names(list_eb_EFPF_fit_NegBinBB) <- paste0("var_fct.", vars_fct_NegBinBB)
  
  list_eb_EFPF_fit_GammaIBP <-  vector(mode = "list", length = length(vars_GammaIBP))
  names(list_eb_EFPF_fit_GammaIBP) <- paste0("var.", vars_GammaIBP)
  
  
  # Loop over the different training set dimensions 
  
  for (j in 1:length(Ns)){
    
    n_train <- Ns[j]
    M <- L - n_train
    
    train_mat <- data_mat[1:n_train,]
    train_mat <- train_mat[, colSums(train_mat) > 0]
    K <- ncol(train_mat)
    
    lab_comb_bb <- paste0("n_train.",n_train,":Nbar.emp")
    lab_comb_ibp <- paste("n_train", n_train, sep = ".")
    
    eb_params_obj_BB <- eb_params(model = "BB", 
                                  init = eb_init_BB, known = eb_known_BB )
    eb_params_obj_IBP <- eb_params(model = "IBP", 
                                  init = eb_init_IBP, known = eb_known_IBP )
    
    # Fit the models
    # PoissonBB
    list_eb_EFPF_fit_PoissonBB[[lab_comb_bb]] <- GibbsFA_eb(feature_matrix = train_mat, 
                                                            model = "PoissonBB", 
                                                            type = "EFPF",
                                                            eb_params =  eb_params_obj_BB)
    
    # NegBinBB
    for (var_fct_NegBinBB in vars_fct_NegBinBB){
      
      list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]][[lab_comb_bb]] <- 
        GibbsFA_eb(feature_matrix = train_mat,
                   model = "NegBinBB", type = "EFPF",
                   eb_params =  eb_params_obj_BB, 
                   var_fct = var_fct_NegBinBB)
      
    }
    
    
    # GammaIBP
    for (var_GammaIBP in vars_GammaIBP){
      
      list_eb_EFPF_fit_GammaIBP[[paste0("var.", var_GammaIBP)]][[lab_comb_ibp]] <-
        GibbsFA_eb(feature_matrix = train_mat,
                   model = "GammaIBP", type = "EFPF",
                   eb_params =  eb_params_obj_IBP,
                   var_GammaIBP = var_GammaIBP)
      
    }
    
    
    # Loop over the different Nbar hyperparameters (only for BB mixtures) 
    
    for (v in 1:length(Nbars)){
      
      Nbar <- Nbars[v]
      lab_comb_bb <- paste0("n_train.",n_train,":Nbar.", Nbar)
      
      # Fit the models
      # PoissonBB
      eb_known_PoissonBB_Nbar <- list("lambda" = Nbar)
      eb_init_PoissonBB_Nbar <- eb_init_BB[! names(eb_init_BB) %in% c("Nhat_prime")]
      eb_params_obj_PoissonBB <- eb_params(model = "PoissonBB",
                                           init = eb_init_PoissonBB_Nbar,
                                           known = eb_known_PoissonBB_Nbar)
      
      
      list_eb_EFPF_fit_PoissonBB[[lab_comb_bb]] <- 
        GibbsFA_eb(feature_matrix = train_mat, 
                   model = "PoissonBB", 
                   type = "EFPF",
                   eb_params =  eb_params_obj_PoissonBB)
      
      # NegBinBB
      for (var_fct_NegBinBB in vars_fct_NegBinBB){
        
        eb_known_NegBinBB_Nbar <- list("mu0" = Nbar, "var_fct" = var_fct_NegBinBB)
        eb_init_NegBinBB_Nbar <- eb_init_BB[! names(eb_init_BB) %in% c("Nhat_prime")]
        eb_params_obj_NegBinBB <- eb_params(model = "NegBinBB",
                                            init = eb_init_NegBinBB_Nbar,
                                            known = eb_known_NegBinBB_Nbar)
        
        list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]][[lab_comb_bb]] <-
          GibbsFA_eb(feature_matrix = train_mat,
                     model = "NegBinBB", type = "EFPF",
                     eb_params =  eb_params_obj_NegBinBB,
                     var_fct = var_fct_NegBinBB)
        
      }
        
      
                                    
    }
    
    
  }
  
  # Save the entire workspace related to the mechanism just performed
  save(list = ls(all.names = TRUE), file =  paste0("R_script_paper/eb_EFPF_",mechanism,"_fit_estimate.RData"))
}




# Main script: EFPF approach -------


# Choose mechanism
mechanism = "beta_pis" # "custom"

# Fit and estimate richness, rarefaction and extrapolation for GibbsFA's (save workspace)
if (!file.exists(paste0("R_script_paper/eb_EFPF_",mechanism,"_fit_estimate.RData"))) {
  
  vars_fct_NegBinBB <- c(10, 1000) 
  vars_GammaIBP <- c(1, 1000)
  
  eb_init_BB <- list(alpha = -1, s = 100, Nhat_prime = 50)
  eb_known_BB <- list()
  
  eb_init_IBP <- list(alpha = 0.5, s = 1, Gamma = 10)
  eb_known_IBP <- list()
  
  # Call the routine to perform simulations
  eb_EFPF_fit_estimate_bounded_scenario(mechanism = mechanism, 
                                        vars_fct_NegBinBB,
                                        vars_GammaIBP,
                                        eb_init_BB, eb_known_BB,
                                        eb_init_IBP, eb_known_IBP,
                                        seed = 123) # seed = 123456 for "custom"
  
}

# Load the Work space
load(paste0("R_script_paper/eb_EFPF_",mechanism,"_fit_estimate.RData"))
Kn <- sapply(Ns, function(n) sum(colSums(data_mat[1:n,]) > 0)  )



## Model-checking on rarefaction ----
n_rare <- Ns[2]
lab_comb_bb <- paste0("n_train.",n_rare,":Nbar.emp")
lab_comb_ibp <- paste0("n_train.",n_rare)

eb_EFPF_fit_PoissonBB_rare <- list_eb_EFPF_fit_PoissonBB[[lab_comb_bb]]
eb_EFPF_fit_NegBinBB_rare <- list_eb_EFPF_fit_NegBinBB[[1]][[lab_comb_bb]]
eb_EFPF_fit_GammaIBP_rare <- list_eb_EFPF_fit_GammaIBP[[1]][[lab_comb_ibp]]

accum_df <- tibble( x = 0:n_rare,
                    n_feat = c(0,rarefaction(data_mat[1:n_rare,], n_reorderings = 20)))

rare_EFPF_PoissonBB <- tibble( lambda_post = unname(unlist(
  rarefaction(object = eb_EFPF_fit_PoissonBB_rare, seed = seed)$lambda_post ))) %>%
  rename(means = lambda_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "PoissonBB",
             x = c(1:n_rare,0))


rare_EFPF_NegBinBB <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_EFPF_fit_NegBinBB_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "NegBinBB",
             x = c(1:n_rare,0))



rare_EFPF_GammaIBP <- tibble( mu0_post = unname(unlist(
  rarefaction(object = eb_EFPF_fit_GammaIBP_rare, seed = seed)$mu0_post ))) %>%
  rename(means = mu0_post) %>%
  add_row(means = 0) %>%
  add_column(Model = "GammaIBP",
             x = c(1:n_rare,0))

# define quantities for plot
rare_EFPF_mixtureBB <- rare_EFPF_NegBinBB %>%
  mutate(Model = "PoissonBB/NegBinBB")

df_rare <- rbind(rare_EFPF_mixtureBB,
                 rare_EFPF_GammaIBP)

df_rare$Model <- factor(df_rare$Model,
                        levels = c( "PoissonBB/NegBinBB", "GammaIBP"))

# for plot
ggplot(accum_df, aes(x = x, y = n_feat)) +
  geom_point(color="black", shape = 18, size = 0.2) + 
  geom_line(data = df_rare, aes(x = x, y = means, color = Model), linetype = "dashed") + 
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) + 
  scale_color_tableau(
    labels = c(
      "PoissonBB/NegBinBB" = "Mixtures of BBs",
      "GammaIBP" = "Mixtures of IBPs"
    )
  ) 


## Model-checking on K_n_r -------
n_knr <- Ns[2]
lab_comb_bb <- paste0("n_train.",n_knr,":Nbar.emp")
lab_comb_ibp <- paste0("n_train.",n_knr)

eb_EFPF_fit_PoissonBB_knr <- list_eb_EFPF_fit_PoissonBB[[lab_comb_bb]]
eb_EFPF_fit_NegBinBB_knr <- list_eb_EFPF_fit_NegBinBB[[1]][[lab_comb_bb]]
eb_EFPF_fit_GammaIBP_knr <- list_eb_EFPF_fit_GammaIBP[[1]][[lab_comb_ibp]]


observed_K_n_r <- tibble( r = 1:n_knr,
                          k_n_r = K_n_r(data_mat[1:n_knr,], n_reorderings = 1)[[paste0('N = ', n_knr)]])


K_n_r_EFPF_PoissonBB <- tibble( lambda_est = unname(unlist(
  K_n_r(object = eb_EFPF_fit_PoissonBB_knr, n = n_knr)[[paste0('N = ', n_knr)]]$lambda_est ))) %>%
  rename(means = lambda_est) %>%
  add_column(Model = "PoissonBB",
             r = 1:n_knr)


K_n_r_EFPF_NegBinBB <- tibble( mu0_est = unname(unlist(
  K_n_r(object = eb_EFPF_fit_NegBinBB_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
  rename(means = mu0_est) %>%
  add_column(Model = "NegBinBB",
             r = 1:n_knr)


K_n_r_EFPF_GammaIBP <- tibble( mu0_est = unname(unlist(
  K_n_r(object = eb_EFPF_fit_GammaIBP_knr, n = n_knr)[[paste0('N = ', n_knr)]]$mu0_est ))) %>%
  rename(means = mu0_est) %>%
  add_column(Model = "GammaIBP",
             r = 1:n_knr)

# define quantities for plot
K_n_r_EFPF_mixtureBB <- K_n_r_EFPF_PoissonBB %>%
  mutate(Model = "PoissonBB/NegBinBB")

df_K_n_r <- rbind(K_n_r_EFPF_mixtureBB,
                  K_n_r_EFPF_GammaIBP)

df_K_n_r$Model <- factor(df_K_n_r$Model, 
                         levels = c("PoissonBB/NegBinBB", "GammaIBP"))

r_positive <- observed_K_n_r %>%
  filter(k_n_r > 0) %>%
  select(r) %>%
  filter(r < 15)

df_K_n_r_plot <- df_K_n_r %>%
  filter(r %in% c(r_positive$r))

observed_K_n_r_plot <- observed_K_n_r %>%
  filter(r %in% c(r_positive$r))

# for plot
ggplot(observed_K_n_r_plot,  aes(x = r, y = k_n_r)) +
  geom_point(color="black", shape = 19, size = 1) +
  geom_line( data = df_K_n_r_plot, aes(x = r, y = means, color = Model), linetype = "dashed") +
  scale_y_log10() +
  xlab("r") + ylab(expression(m[r])) + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) + 
  scale_color_tableau(
    labels = c(
    "PoissonBB/NegBinBB" = "Mixtures of BBs",
    "GammaIBP" = "Mixtures of IBPs"
    )
  )



## Prediction: richness and extrapolation -----------

### Richness -------

# 1) Plot Richness: Point plot expected value (number of features)
labels_comb_bb <- paste(rep(paste("n_train", Ns, sep = "."), each = length(Nbars)+1),
                        c("Nbar.emp" , paste("Nbar", Nbars, sep = ".")), sep=":")

# PoissonBB
richness_EFPF_PoissonBB_df <- tibble(estimate = unname(sapply(list_eb_EFPF_fit_PoissonBB, function(x)
  total_richness(x)$lambda_post + ncol(x$feature_matrix)) ) ) %>%
  add_column(Model = "Poisson BB", 
             Nbar = rep(c("EB", Nbars), length(Ns)),
             n_train = rep(Ns, each = length(Nbars)+ 1),
             n_train_idx = rep(c(1,2,3), each = length(Nbars)+ 1)) 

# NegBinBB
richness_EFPF_NegBinBB_df <- tibble(estimate = numeric(), Model = character(),
                                    Nbar = character(), 
                                    n_train = integer(), n_train_idx = integer())

for (var_fct_NegBinBB in vars_fct_NegBinBB){
  list_eb_EFPF_fit_NegBinBB_var <- list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]]
  
  richness_EFPF_NegBinBB_df_var <- tibble(estimate = unname(sapply(list_eb_EFPF_fit_NegBinBB_var, function(x)
    total_richness(x)$mu0_post + ncol(x$feature_matrix)) ) ) %>%
    add_column(Model = paste0("NegBinomial BB x",var_fct_NegBinBB) , 
               Nbar = rep(c("EB", Nbars), length(Ns)),
               n_train = rep(Ns, each = length(Nbars)+ 1),
               n_train_idx = rep(c(1,2,3), each = length(Nbars)+ 1))
  
  richness_EFPF_NegBinBB_df <- bind_rows(richness_EFPF_NegBinBB_df, richness_EFPF_NegBinBB_df_var)
  
}


joint_richness_long <- bind_rows(richness_EFPF_PoissonBB_df,
                                 richness_EFPF_NegBinBB_df) %>%
  mutate(Model = fct_relevel(Model, c("Poisson BB", 
                             paste0("NegBinomial BB x", vars_fct_NegBinBB))) )


joint_richness_long <- joint_richness_long %>%
  mutate(n_train_latex = paste("Scenario~", LETTERS[n_train_idx], ":", "~n == ", n_train, "~','~K[n] == ", Kn[n_train_idx], sep = ""))


# for plot
ggplot(joint_richness_long, aes( y=estimate, x=Model, shape = Nbar)) +
  geom_point(size = 2) +
  facet_wrap(~ n_train_latex,
             labeller = label_parsed,
             scales = "free_x", nrow = 1) +
  theme_light() +
  geom_hline(aes(yintercept = H), linetype = "dashed") +
  theme(legend.position = "top") +
  ylab("Posterior mean of N") +
  scale_y_continuous(breaks = pretty_breaks()) +
  rremove("xlab") +
  scale_shape_discrete(name = "Prior mean of N") +
  theme(aspect.ratio = 1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# 2) Plot Richness: posterior distributions for given Nbar  
Nbar_plot <- 400

# PoissonBB
list_eb_EFPF_est_PoissonBB <- list_eb_EFPF_fit_PoissonBB[grepl("Nbar.emp", names(list_eb_EFPF_fit_PoissonBB) )]
params_richness_EFPF_est_PoissonBB <- tibble( lambda_prime = unname(sapply(list_eb_EFPF_est_PoissonBB, function(x)
  total_richness(x)$lambda_post))) %>%
  add_column(n_train = Ns,
             n_train_idx = rep(c(1,2,3), each = 1),
             Model = "Poisson BB", PM = "est") %>%
  mutate(lb = qpois(0.025, lambda_prime, lower.tail = TRUE, log.p = FALSE),
         ub = qpois(0.975, lambda_prime, lower.tail = TRUE, log.p = FALSE) )

list_eb_EFPF_fixed_PoissonBB <- list_eb_EFPF_fit_PoissonBB[grepl(paste0("Nbar.",Nbar_plot), names(list_eb_EFPF_fit_PoissonBB) )]
params_richness_EFPF_fixed_PoissonBB <- tibble( lambda_prime = unname(sapply(list_eb_EFPF_fixed_PoissonBB, function(x)
  total_richness(x)$lambda_post))) %>%
  add_column(n_train = Ns,
             n_train_idx = rep(c(1,2,3), each = 1),
             Model = "Poisson BB", PM = "fixed") %>%
  mutate(lb = qpois(0.025, lambda_prime, lower.tail = TRUE, log.p = FALSE),
         ub = qpois(0.975, lambda_prime, lower.tail = TRUE, log.p = FALSE) )

# NegBinBB
params_richness_EFPF_est_NegBinBB <- tibble( n0_prime = numeric(),
                                             mu0_prime = numeric(),
                                             n_train = integer(),
                                             n_train_idx = integer(),
                                             Model = character(),
                                             PM = character(),
                                             p_prime = numeric(),
                                             lb = numeric(), ub = numeric())
params_richness_EFPF_fixed_NegBinBB <- tibble( n0_prime = numeric(),
                                             mu0_prime = numeric(),
                                             n_train = integer(),
                                             n_train_idx = integer(),
                                             Model = character(),
                                             PM = character(),
                                             p_prime = numeric(),
                                             lb = numeric(), ub = numeric())

for (var_fct_NegBinBB in vars_fct_NegBinBB){
  
  # est
  list_eb_EFPF_est_NegBinBB_var <- 
    list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]][grepl("Nbar.emp", names(list_eb_EFPF_fit_NegBinBB[[1]]) )]
  
  params_richness_EFPF_est_NegBinBB_var <- tibble( n0_prime = unname(sapply(list_eb_EFPF_est_NegBinBB_var, function(x)
    total_richness(x)$n0_post)),
    mu0_prime = unname(sapply(list_eb_EFPF_est_NegBinBB_var, function(x)
      total_richness(x)$mu0_post))) %>%
    add_column(n_train = Ns,
               n_train_idx = rep(c(1,2,3), each = 1),
               Model = paste0("NegBinomial BB x", var_fct_NegBinBB),
               PM = "est") %>%
    mutate(p_prime = 1/(mu0_prime/n0_prime + 1),
           lb = qnbinom(0.025, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE),
           ub = qnbinom(0.975, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE) )
  
  params_richness_EFPF_est_NegBinBB <- bind_rows(params_richness_EFPF_est_NegBinBB,
                                                 params_richness_EFPF_est_NegBinBB_var)

  
  #fixed
  list_eb_EFPF_fixed_NegBinBB_var <- 
    list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]][grepl(paste0("Nbar.",Nbar_plot), names(list_eb_EFPF_fit_NegBinBB[[1]]) )]
  
  params_richness_EFPF_fixed_NegBinBB_var <- tibble( n0_prime = unname(sapply(list_eb_EFPF_fixed_NegBinBB_var, function(x)
    total_richness(x)$n0_post)),
    mu0_prime = unname(sapply(list_eb_EFPF_fixed_NegBinBB_var, function(x)
      total_richness(x)$mu0_post))) %>%
    add_column(n_train = Ns,
               n_train_idx = rep(c(1,2,3), each = 1),
               Model = paste0("NegBinomial BB x", var_fct_NegBinBB),
               PM = "fixed") %>%
    mutate(p_prime = 1/(mu0_prime/n0_prime + 1),
           lb = qnbinom(0.025, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE),
           ub = qnbinom(0.975, size = n0_prime, prob = p_prime, lower.tail = TRUE, log.p = FALSE) )
  
  
  params_richness_EFPF_fixed_NegBinBB <- bind_rows(params_richness_EFPF_fixed_NegBinBB,
                                                   params_richness_EFPF_fixed_NegBinBB_var)

  
  
}



bounds_distr <- tibble( lb = min(params_richness_EFPF_est_PoissonBB$lb + Kn, 
                                 params_richness_EFPF_est_NegBinBB$lb + Kn,
                                 params_richness_EFPF_fixed_PoissonBB$lb + Kn, 
                                 params_richness_EFPF_fixed_NegBinBB$lb + Kn,
                                 H),
                        ub = max(params_richness_EFPF_est_PoissonBB$ub + Kn,
                                 params_richness_EFPF_est_NegBinBB$ub + Kn,
                                 params_richness_EFPF_fixed_PoissonBB$ub + Kn,
                                 params_richness_EFPF_fixed_NegBinBB$ub + Kn,
                                 H))

# density computation
# PoissonBB
dens_richness_EFPF_est_PoissonBB <- bind_rows(lapply(1:length(Ns), function(j) 
  tibble( x = bounds_distr$lb: bounds_distr$ub) %>%
    mutate( y = dpois(x - Kn[j] , lambda = params_richness_EFPF_est_PoissonBB$lambda_prime[j])) %>%
    add_column(Model = "Poisson BB", PM = "est",
               n_train = params_richness_EFPF_est_PoissonBB$n_train[j],
               n_train_idx = params_richness_EFPF_est_PoissonBB$n_train_idx[j])
))

dens_richness_EFPF_fixed_PoissonBB <- bind_rows(lapply(1:length(Ns), function(j) 
  tibble( x = bounds_distr$lb: bounds_distr$ub) %>%
    mutate( y = dpois(x - Kn[j] , lambda = params_richness_EFPF_fixed_PoissonBB$lambda_prime[j])) %>%
    add_column(Model = "Poisson BB", PM = "fixed",
               n_train = params_richness_EFPF_fixed_PoissonBB$n_train[j],
               n_train_idx = params_richness_EFPF_fixed_PoissonBB$n_train_idx[j] )
))

# NegBinBB
dens_richness_EFPF_est_NegBinBB <- tibble(x = integer(), y = numeric(),
                                          Model = character(),
                                          PM = character(),
                                          n_train = integer(),
                                          n_train_idx = integer())

dens_richness_EFPF_fixed_NegBinBB <- tibble(x = integer(), y = numeric(),
                                            Model = character(),
                                            PM = character(),
                                            n_train = integer(),
                                            n_train_idx = integer())


for (var_fct_NegBinBB in vars_fct_NegBinBB){
  
  # est
  params_richness_EFPF_est_NegBinBB_var <- params_richness_EFPF_est_NegBinBB %>%
    filter(Model == paste0("NegBinomial BB x", var_fct_NegBinBB))
  
  dens_richness_EFPF_est_NegBinBB_var <- bind_rows(lapply(1:length(Ns), function(j) 
    tibble( x = bounds_distr$lb: bounds_distr$ub) %>%
      mutate( y = dnbinom(x - Kn[j], size = params_richness_EFPF_est_NegBinBB_var$n0_prime[j], 
                          prob = params_richness_EFPF_est_NegBinBB_var$p_prime[j])) %>%
      add_column(Model = paste0("NegBinomial BB x", var_fct_NegBinBB), 
                 PM = "est",
                 n_train = params_richness_EFPF_est_NegBinBB_var$n_train[j],
                 n_train_idx = params_richness_EFPF_est_NegBinBB_var$n_train_idx[j])
  ))
  
  dens_richness_EFPF_est_NegBinBB <- bind_rows(dens_richness_EFPF_est_NegBinBB,
                                               dens_richness_EFPF_est_NegBinBB_var)
  
  
  # fixed
  params_richness_EFPF_fixed_NegBinBB_var <- params_richness_EFPF_fixed_NegBinBB %>%
    filter(Model == paste0("NegBinomial BB x", var_fct_NegBinBB))
  
  dens_richness_EFPF_fixed_NegBinBB_var <- bind_rows(lapply(1:length(Ns), function(j) 
    tibble( x = bounds_distr$lb: bounds_distr$ub) %>%
      mutate( y = dnbinom(x - Kn[j], size = params_richness_EFPF_fixed_NegBinBB_var$n0_prime[j], 
                          prob = params_richness_EFPF_fixed_NegBinBB_var$p_prime[j])) %>%
      add_column(Model = paste0("NegBinomial BB x", var_fct_NegBinBB),
                 PM = "fixed",
                 n_train = params_richness_EFPF_fixed_NegBinBB_var$n_train[j],
                 n_train_idx = params_richness_EFPF_fixed_NegBinBB_var$n_train_idx[j])
  ))
  
  dens_richness_EFPF_fixed_NegBinBB <- bind_rows(dens_richness_EFPF_fixed_NegBinBB,
                                                 dens_richness_EFPF_fixed_NegBinBB_var)
  
}


# Set Ns and Kn to plot
idx_plot <- c(1,2,3) 
Ns_plot <- Ns[idx_plot] 
Kn_plot <- Kn[idx_plot]

dens_richnesses <- rbind(dens_richness_EFPF_est_PoissonBB,
                         dens_richness_EFPF_est_NegBinBB,
                         dens_richness_EFPF_fixed_PoissonBB,
                         dens_richness_EFPF_fixed_NegBinBB ) %>%
  filter(n_train %in% Ns_plot)


dens_richnesses$Model <- factor(dens_richnesses$Model, 
                                levels = c("Poisson BB",
                                           paste0("NegBinomial BB x", vars_fct_NegBinBB)))

dens_richnesses$PM <- factor(dens_richnesses$PM)
levels(dens_richnesses$PM) <- c("EB", "Bayesian")
dens_richnesses <- dens_richnesses %>%
  mutate(n_train_latex = paste("Scenario~", LETTERS[n_train_idx], ":", "~n == ", n_train, "~','~K[n] == ", Kn[n_train_idx], sep = ""))

# for plot
ggplot(dens_richnesses, aes(x = x, y = y, color = Model)) +
  geom_line() +
  geom_vline(aes(xintercept = H), linetype="dashed") +
  facet_grid(PM ~ n_train_latex,
             labeller = label_parsed,
             scales = "free") +
  theme_light() +
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  xlab("# distinct features") + ylab("Probability") + 
  scale_color_tableau()



### Extrapolation -----

# read or compute GT and Chao estimators
if (!file.exists(paste0("R_script_paper/eb_Freq_",mechanism,"_fit_estimate.RData"))) {
  
  M <- 300 
  
  # List to store Chao's estimates
  list_rare_extr_Chao <- vector(mode="list")
  # List to store smoothed GT's estimates
  list_extr_GT <- vector(mode="list")
  
  for (j in 1:length(Ns)){
    
    n_train <- Ns[j]
    
    train_mat <- data_mat[1:n_train,]
    train_mat <- train_mat[, colSums(train_mat) > 0]
    K <- ncol(train_mat)
    
    lab_comb_bb <- paste0("n_train.",n_train,":Nbar.emp")
    lab_comb_ibp <- paste("n_train", n_train, sep = ".")
    
    # A) Chao's rarefaction and extrapolation
    
    if (j < 3){
      # Determine the frequency vector of the training sets
      Q_vec <- colSums(train_mat)
      Q_vec <- Q_vec[Q_vec>0]
      
      # Compute the curves with confidence intervals
      fit_Chao <- iNEXT.Sam(Spec = Q_vec, T = n_train, endpoint = n_train + M)
      
      rare_Chao <- as_tibble(fit_Chao[["q=0"]]) %>%
        select(-Cov.hat) %>%
        rename(medians = D0.hat, lbs = Norm.CI.Low, ubs = Norm.CI.High)
      
      list_rare_extr_Chao[[lab_comb_ibp]] <- as.data.frame(rare_Chao)
      
    }
    
    
    # B) Smoothed Good-Toulmin extrapolation
    
    # Compute SFS vector and CTS vector
    sfs <- tabulate(colSums(train_mat))
    cts <- sapply(2:n_train, function(n) ncol(train_mat[1:n,colSums(train_mat[1:n,]) > 0])   )
    cts <- c(0, sum(train_mat[1,]) , cts)
    
    list_extr_GT[[lab_comb_ibp]] <- predict_good_toulmin(n_train, M, sfs, cts, alternative = 0)$preds
    
    
  }
  
  # Save the entire workspace related to the mechanism just performed
  save(list = ls(all.names = TRUE), file =  paste0("R_script_paper/eb_Freq_",mechanism,"_fit_estimate.RData"))
  
}

# Load the Work space
load(paste0("R_script_paper/eb_Freq_",mechanism,"_fit_estimate.RData"))



# Extract accumulation curve of the observed sample (or average accumulation)

accum_df <- tibble(n_feat = unlist(sapply(Ns, function(n) 
  c(0,rarefaction(data_mat[1:(n + M), ], n_reorderings = 1)))),
  n_train = rep(Ns, times = Ns + M +1),
  n_train_idx = rep(c(1,2,3), times = Ns + M +1 ),
  type = unlist(sapply(Ns, function(n) c(rep("train", n +1), rep("test", M)))),
  x = unlist(sapply(Ns, function(n) 0:(n + M))))

accum_df_train <- accum_df %>%
  filter(type == "train")
accum_df_test <- accum_df %>%
  filter(type == "test")


# Poisson
list_eb_EFPF_PoissonBB <- list_eb_EFPF_fit_PoissonBB[grepl("Nbar.emp", names(list_eb_EFPF_fit_PoissonBB) )]
extr_EFPF_PoissonBB_df <- tibble(lambda = unname(unlist(lapply(list_eb_EFPF_PoissonBB, function(x)
  extrapolation(object = x, M = M, seed = seed)$lambda_post ))),
  n_train = rep(Ns, each = M),
  n_train_idx = rep(c(1,2,3), each = M),
  Kn = rep(Kn, each = M)) %>%
  mutate(lb = qpois(0.025, lambda, lower.tail = TRUE, log.p = FALSE),
         ub = qpois(0.975, lambda, lower.tail = TRUE, log.p = FALSE)) %>%
  rename(means = lambda) %>%
  add_row(means = 0, lb = 0, ub = 0, n_train = Ns, n_train_idx = c(1,2,3), Kn = Kn) %>%
  mutate(means = means + Kn, lb = lb + Kn, ub = ub + Kn) %>%
  add_column(x = c(unlist(sapply(Ns, function(n) (n+1):(M+n))), Ns),
             Model = "Poisson BB")

extr_EFPF_PoissonBB_df$x <- as.integer(extr_EFPF_PoissonBB_df$x)
extr_EFPF_PoissonBB_df <- extr_EFPF_PoissonBB_df %>%
  select(means, lb, ub, n_train, n_train_idx, x, Model)


# NegBin
extr_EFPF_NegBinBB_df <- tibble(means = numeric(), 
                                lb = numeric(), ub = numeric(),
                                n_train = integer(), 
                                n_train_idx = integer(), 
                                x = integer(), Model = character())

for (var_fct_NegBinBB in vars_fct_NegBinBB){
  
  list_eb_EFPF_NegBinBB_var <- 
    list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]][grepl("Nbar.emp", names(list_eb_EFPF_fit_NegBinBB[[1]]) )]
  
  extr_EFPF_NegBinBB_df_var <- tibble(mu0 = unname(unlist(lapply(list_eb_EFPF_NegBinBB_var, function(x)
    extrapolation(object = x, M = M, seed = seed)$mu0_post ))),
    n0 = unname(unlist(lapply(list_eb_EFPF_NegBinBB_var, function(x)
      extrapolation(object = x, M = M, seed = seed)$n0_post ))),
    n_train = rep(Ns, each = M),
    n_train_idx = rep(c(1,2,3), each = M),
    Kn = rep(Kn, each = M)) %>%
    mutate(p = 1/(mu0/n0 + 1),
           lb = qnbinom(0.025, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE),
           ub = qnbinom(0.975, size = n0, prob = p, lower.tail = TRUE, log.p = FALSE)) %>%
    rename(means = mu0) %>%
    add_row(means = 0, lb = 0, ub = 0, n_train = Ns, n_train_idx = c(1,2,3), Kn = Kn) %>%
    mutate(means = means + Kn, lb = lb + Kn, ub = ub + Kn) %>%
    add_column(x = c(unlist(sapply(Ns, function(n) (n+1):(M+n))), Ns),
               Model = paste0("NegBinomial BB x", var_fct_NegBinBB))
  
  extr_EFPF_NegBinBB_df_var$x <- as.integer(extr_EFPF_NegBinBB_df_var$x)
  extr_EFPF_NegBinBB_df_var <- extr_EFPF_NegBinBB_df_var %>%
    select(means, lb, ub, n_train, n_train_idx, x, Model)
  
  extr_EFPF_NegBinBB_df <- bind_rows(extr_EFPF_NegBinBB_df, 
                                     extr_EFPF_NegBinBB_df_var)
  
}



# Set Ns and Kn to plot
idx_plot <- c(1,2) 
Ns_plot <- Ns[idx_plot] 
Kn_plot <- Kn[idx_plot]

df_extr_GT_long <- list_extr_competitor_to_long(list_extr_GT, model = "GT") %>%
  rename(Model = model) %>%
  add_column(n_train_idx = rep(c(1,2,3), each = M + 1)) %>%
  filter(t < n_train + M + 1, n_train %in% Ns_plot)
df_extr_Chao_long <- list_extr_competitor_to_long(list_rare_extr_Chao, model = "Chao") %>%
  rename(Model = model) %>%
  filter(t < n_train + M + 1, n_train %in% Ns_plot) %>%
  mutate(n_train_idx = match(n_train, Ns_plot))

temp <- tibble(n_train_latex = paste("Scenario~", LETTERS[c(1,2,3)], ":", "~n == ", Ns, "~','~K[n] == ", Kn, sep = ""), 
               xvalues = Ns) %>%
  filter(xvalues %in% Ns_plot)

extr_EFPF_PoissonBB_df <- extr_EFPF_PoissonBB_df %>%
  add_column(Model_gen = "PoissonBB/NegBinBB")
extr_EFPF_NegBinBB_df <- extr_EFPF_NegBinBB_df %>%
  add_column(Model_gen = "PoissonBB/NegBinBB")

extr_all_df <- rbind(extr_EFPF_PoissonBB_df, 
                     extr_EFPF_NegBinBB_df) %>%
  filter(n_train %in% Ns_plot)
                     

extr_all_df$Model <- factor(extr_all_df$Model, 
                            levels = c("Poisson BB", 
                                       paste0("NegBinomial BB x", vars_fct_NegBinBB),
                                       paste0("GammaIBP, var = ", vars_GammaIBP)))

accum_df_train <- accum_df_train %>%
  filter(n_train %in% Ns_plot)
accum_df_test <- accum_df_test %>%
  filter(n_train %in% Ns_plot)


accum_df_train <- accum_df_train %>%
  mutate(n_train_latex = paste("Scenario~", LETTERS[n_train_idx], ":", "~n == ", n_train, "~','~K[n] == ", Kn[n_train_idx], sep = ""))

accum_df_test <- accum_df_test %>%
  mutate(n_train_latex = paste("Scenario~", LETTERS[n_train_idx], ":", "~n == ", n_train, "~','~K[n] == ", Kn[n_train_idx], sep = ""))

extr_all_df <- extr_all_df %>%
  mutate(n_train_latex = paste("Scenario~", LETTERS[n_train_idx], ":", "~n == ", n_train, "~','~K[n] == ", Kn[n_train_idx], sep = ""))

df_extr_GT_long <- df_extr_GT_long %>%
  mutate(n_train_latex = paste("Scenario~", LETTERS[n_train_idx], ":", "~n == ", n_train, "~','~K[n] == ", Kn[n_train_idx], sep = ""))

df_extr_Chao_long <- df_extr_Chao_long %>%
  mutate(n_train_latex = paste("Scenario~", LETTERS[n_train_idx], ":", "~n == ", n_train, "~','~K[n] == ", Kn[n_train_idx], sep = ""))

# for plot
plot_means <- ggplot() +
  geom_point(data = accum_df_train, aes(x = x, y = n_feat),
             color="black", shape = 19, size = 0.1) +
  facet_grid(~ n_train_latex,
             labeller = label_parsed,
             scales = "free_x")  +
  geom_point( data = accum_df_test, aes(x = x, y = n_feat),
              color="black", shape = 19, size = 0.1) +
  geom_line(data = extr_all_df, aes(x = x, y = means, color = Model),
            linetype = "dashed") +
  geom_line(data = df_extr_GT_long, aes(t, value, color = Model), linetype = "dashed") +
  geom_line(data = df_extr_Chao_long, aes(t, medians, color = Model), linetype = "dashed") +
  geom_vline(data = temp, mapping =  aes(xintercept = xvalues) , linetype = "dashed", color = "grey") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) +
  scale_color_tableau()

plot_ribbons <- ggplot() +
  geom_point(data = accum_df_train, aes(x = x, y = n_feat),
             color="black", shape = 19, size = 0.1) +
  geom_point( data = accum_df_test, aes(x = x, y = n_feat),
              color="black", shape = 19, size = 0.1) +
  geom_ribbon(data = extr_all_df, aes(x = x, ymin = lb, ymax = ub, fill = Model), color = NA, alpha = 0.4) +
  scale_fill_manual(values = c("Poisson BB" = "grey10", "NegBinomial BB x10" = "grey50", "NegBinomial BB x1000" = "grey80")) +
  facet_grid(~ n_train_latex,
             labeller = label_parsed,
             scales = "free_x")  +
  geom_vline(data = temp, mapping =  aes(xintercept = xvalues) , linetype = "dashed", color = "grey") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) +
  scale_color_tableau()
  

combined_plot <- plot_means / plot_ribbons
combined_plot




