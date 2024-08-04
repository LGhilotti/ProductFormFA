#
# Unbounded scenario simulations ####
#

rm(list=ls())
library(ProductFormFA)
library(tidyverse)
library(scales)
library(ggpubr)
library(ggthemes)
library(patchwork)

source("R_script_paper/Routine_Chao.R")
source("R_script_paper/utils.R")


# Function to generate binary matrix for scenario B ----

generate_data <- function(xi, n, H, seed = 1234){
  
  set.seed(seed)
  
  pis <- 1/((1:H)^xi)
  
  data_mat <- matrix(rbinom(n*H, size = 1, prob = rep(pis, n)),
                     nrow = n, ncol = H, byrow = T )
  
  return(data_mat)
  
}


# Function to produce fit and estimate on specific polynomial scenario ------

eb_EFPF_fit_estimate_polynomial_scenario <- function(xi, 
                                                     vars_fct_NegBinBB,
                                                     vars_GammaIBP,
                                                     eb_init_BB, eb_known_BB,
                                                     eb_init_IBP, eb_known_IBP,
                                                     seed = 1234){
  
  if (xi <= 0){
    stop("Invalid generating mechanism.")
  }
  
  # Set training dimensions and dimension of the whole sample 
  Ns <- c(10, 50, 250) 
  N_max <- Ns[length(Ns)]
  L <- N_max + 400
  
  # Set maximum number of features
  H <- 10^6
  
  # Generate data 
  data_mat <- generate_data(xi = xi, n = L, H = H, seed = seed)
  data_mat <- data_mat[, colSums(data_mat) > 0]
  
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
    
  }
  
  # Save the entire workspace related to the mechanism just performed
  save(list = ls(all.names = TRUE), file =  paste0("R_script_paper/eb_EFPF_poly_",xi,"_fit_estimate.RData"))
}



# Main script: EFPF approach -------


# Choose mechanism
xi = 1 

# Fit and estimate richness, rarefaction and extrapolation for GibbsFA's (save workspace)
if (!file.exists(paste0("R_script_paper/eb_EFPF_poly_",xi,"_fit_estimate.RData"))) {
  
  vars_fct_NegBinBB <- c(10, 1000) 
  vars_GammaIBP <- c(0.01, 100)
  
  eb_init_BB <- list(alpha = -10, s = 10, Nhat_prime = 100)
  eb_known_BB <- list()
  
  eb_init_IBP <- list(alpha = 0.5, s = 1, Gamma = 10)
  eb_known_IBP <- list()
  
  # Call the routine to perform simulations
  eb_EFPF_fit_estimate_polynomial_scenario(xi = xi, 
                                           vars_fct_NegBinBB,
                                           vars_GammaIBP,
                                           eb_init_BB, eb_known_BB,
                                           eb_init_IBP, eb_known_IBP,
                                           seed = 123456)
  
}

# Load the Work space
load(paste0("R_script_paper/eb_EFPF_poly_",xi,"_fit_estimate.RData"))
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
                        levels = c("PoissonBB/NegBinBB", "GammaIBP"))

# for plot
ggplot(accum_df, aes(x = x, y = n_feat)) +
  geom_point(color="black", shape = 19, size = 0.4) + 
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
  filter(r < 10)

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



## Extrapolation ----

# read or compute GT estimator
if (!file.exists(paste0("R_script_paper/eb_Freq_poly_",xi,"_fit_estimate.RData"))) {
  
  M <- 400
  
  # List to store smoothed GT's estimates
  list_extr_GT <- vector(mode="list")
  
  for (j in 1:length(Ns)){
    
    n_train <- Ns[j]
    
    train_mat <- data_mat[1:n_train,]
    train_mat <- train_mat[, colSums(train_mat) > 0]
    K <- ncol(train_mat)
    
    lab_comb_bb <- paste0("n_train.",n_train,":Nbar.emp")
    lab_comb_ibp <- paste("n_train", n_train, sep = ".")
    
    # Smoothed Good-Toulmin extrapolation
    
    # Compute SFS vector and CTS vector
    sfs <- tabulate(colSums(train_mat))
    cts <- sapply(2:n_train, function(n) ncol(train_mat[1:n,colSums(train_mat[1:n,]) > 0])   )
    cts <- c(0, sum(train_mat[1,]) , cts)
    
    list_extr_GT[[lab_comb_ibp]] <- predict_good_toulmin(n_train, M, sfs, cts, alternative = 0)$preds
    
    
  }
  
  # Save the entire workspace related to the mechanism just performed
  save(list = ls(all.names = TRUE), file =  paste0("R_script_paper/eb_Freq_poly_",xi,"_fit_estimate.RData"))
  
}

# Load the Work space
load(paste0("R_script_paper/eb_Freq_poly_",xi,"_fit_estimate.RData"))


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


# GammaIBP
extr_EFPF_GammaIBP_df <- tibble(means = numeric(), 
                                lb = numeric(), ub = numeric(),
                                n_train = integer(), 
                                n_train_idx = integer(),
                                x = integer(), Model = character())

for (var_GammaIBP in vars_GammaIBP){
  
  list_eb_EFPF_GammaIBP_var <- 
    list_eb_EFPF_fit_GammaIBP[[paste0("var.", var_GammaIBP)]]
  
  extr_EFPF_GammaIBP_df_var <- tibble(mu0 = unname(unlist(lapply(list_eb_EFPF_GammaIBP_var, function(x)
    extrapolation(object = x, M = M, seed = seed)$mu0_post ))),
    n0 = unname(unlist(lapply(list_eb_EFPF_GammaIBP_var, function(x)
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
               Model = paste0("Gamma IBP, Variance: ", var_GammaIBP))
  
  extr_EFPF_GammaIBP_df_var$x <- as.integer(extr_EFPF_GammaIBP_df_var$x)
  extr_EFPF_GammaIBP_df_var <- extr_EFPF_GammaIBP_df_var %>%
    select(means, lb, ub, n_train, n_train_idx, x, Model)
  
  
  extr_EFPF_GammaIBP_df <- bind_rows(extr_EFPF_GammaIBP_df, 
                                     extr_EFPF_GammaIBP_df_var)
  
}

# Set Ns and Kn to plot
idx_plot <- c(1,2,3) 
Ns_plot <- Ns[idx_plot] 
Kn_plot <- Kn[idx_plot]

df_extr_GT_long <- list_extr_competitor_to_long(list_extr_GT, model = "GT") %>%
  rename(Model = model) %>%
  add_column(n_train_idx = rep(c(1,2,3), each = M + 1)) %>%
  filter(t < n_train + M + 1, n_train %in% Ns_plot)


temp <- tibble(n_train_latex = paste("Scenario~", LETTERS[c(1,2,3)], ":", "~n == ", Ns, "~','~K[n] == ", Kn, sep = ""), 
               xvalues = Ns) %>%
  filter(xvalues %in% Ns_plot)


extr_all_df <- rbind(extr_EFPF_GammaIBP_df)  %>%
  filter(n_train %in% Ns_plot)

extr_all_df$Model <- factor(extr_all_df$Model, 
                            levels = paste0("Gamma IBP, Variance: ", vars_GammaIBP))
                                       
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
  scale_fill_manual(values = c(
    "Gamma IBP, Variance: 0.01" = "grey10",
    "Gamma IBP, Variance: 100" = "grey60")) +
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

