############ 1) Generate the data #####################

# set number of individuals
L <- 1500 # 1500
N <- 50 # dimension of the training (used as whole sample when we don't consider test/train)
M <- L - N
# set maximum number of features in the population
H <- 200

# value of xi
xis = c(0.6, 0.8, 1, 1.2, 1.4)

seed = 1234

gg_list_poiss <- vector(mode="list", length = length(xis))
gg_list_negbin <- vector(mode="list", length = length(xis))

for (j in 1:length(xis)){
  xi <- xis[j] 
  
  # generate data from zipf, in form of binary matrix
  data_mat <- rzipf(n = L, K = H, xi = xi, seed = seed)
  
  # convert the binary matrix into list of features
  data_list <- create_features_list(data_mat)
  plot_trajectory(data_list)
  
  set.seed(1234)
  idx_train <- sample(1:L, N)
  train_list <- data_list[idx_train]
  test_list <- data_list[setdiff(1:L, idx_train)]
  plot_trajectory(train_list)
  
  train_mat <- create_features_matrix(train_list)
  test_mat <- create_features_matrix(test_list)
  
  
  ############ 2) Set parameters for the 3 models ###############
  
  ########### 2.1) Set parameters for BB (Poisson and Neg-Bin)
  
  # Set prior hyperparameters specific for poisson
  a_l_poiss <- 50
  b_l_poiss <- 0.1
  print(paste0("E(lambda) = ", a_l_poiss/ b_l_poiss))
  print(paste0("Var(lambda) = ", a_l_poiss/ (b_l_poiss^2)))
  
  # Set prior hyperparameters specific for NB
  q_star_nb <- 0.01
  print(paste0("E(nstar) = ", 1/ q_star_nb))
  print(paste0("Var(nstar) = ", (1-q_star_nb)/ (q_star_nb^2)) )
  alpha_p_nb <- 0.1
  beta_p_nb <- 0.1
  print(paste0("E(p) = ", alpha_p_nb/(alpha_p_nb + beta_p_nb)))
  print(paste0("Var(p) = ",(alpha_p_nb*beta_p_nb)/((alpha_p_nb + beta_p_nb)^2 *(alpha_p_nb+beta_p_nb+1)) ))
  
  # Set other hyperparameters
  a_alpha_bb <- 10
  b_alpha_bb <- 1
  print(paste0("E(alpha_bar) = ", a_alpha_bb/ b_alpha_bb))
  print(paste0("Var(alpha_bar) = ", a_alpha_bb/ (b_alpha_bb^2)))
  a_s_bb <- 2
  b_s_bb <- 0.2
  print(paste0("E(s) = ", a_s_bb/ b_s_bb))
  print(paste0("Var(s) = ", a_s_bb/ (b_s_bb^2)))
  
  # Set initial values for poisson 
  lambda_0_poiss <- 100
  # Set initial values for NB
  nstar_0_nb <- 100
  p_0_nb <- 0.2
  # Set initial values for other parameters
  alpha_bar_0_bb <- 1
  s_0_bb <- 1
  
  
  ########## 2.2) Set parameters for the Gamma IBP
  
  # Set prior hyperparameters
  p_ibp <- 0.05
  print(paste0("E(a) = ", 1/ p_ibp))
  print(paste0("Var(a) = ", (1-p_ibp)/ (p_ibp^2)) )
  r_ibp <- 1
  t_ibp <- 0.1
  print(paste0("E(b) = ", r_ibp/ t_ibp))
  print(paste0("Var(b) = ", r_ibp/ (t_ibp^2)))
  a_alpha_ibp <- 0.5
  b_alpha_ibp <- 0.5
  print(paste0("E(alpha) = ", a_alpha_ibp/ (a_alpha_ibp + b_alpha_ibp)))
  print(paste0("Var(alpha) = ", a_alpha_ibp*b_alpha_ibp/ (a_alpha_ibp + b_alpha_ibp)^2 /(a_alpha_ibp+b_alpha_ibp+1)))
  a_s_ibp <- 2
  b_s_ibp <- 0.2
  print(paste0("E(s) = ", a_s_ibp/ b_s_ibp))
  print(paste0("Var(s) = ", a_s_ibp/ (b_s_ibp^2)))
  
  print(paste0("E(theta) = ", a_s_ibp/ b_s_ibp -  a_alpha_ibp/ (a_alpha_ibp + b_alpha_ibp)))
  print(paste0("Var(theta) = ", a_s_ibp/ (b_s_ibp^2) + 
                 a_alpha_ibp*b_alpha_ibp/ (a_alpha_ibp + b_alpha_ibp)^2 /(a_alpha_ibp+b_alpha_ibp+1)))
  
  
  # Set initial values for the parameters of Gamma IBP
  a_0_ibp <- 5
  b_0_ibp <- 1
  alpha_0_ibp <- 0.5
  s_0_ibp <- 15
  
  
  ########## 2.3) Set MCMC parameters (common to all 3 models)
  
  S_poiss <- S_negbin <- 4*10^4
  S_ibp <- 4*10^4
  n_burnin_poiss <- n_burnin_negbin <- 5*10^3
  n_burnin_ibp <- 5*10^3
  thin_poiss <- thin_negbin <- thin_ibp <- 5
  seed <- 1234
  
  
  ################# 3) Run BB + Poisson ###########
  
  # Set tau for MALA
  tau_poiss <- 0.0005
  
  output_poiss <- gibbs_sampler_poiss(Z = train_mat,
                                      lambda_0_poiss, alpha_bar_0_bb, s_0_bb,
                                      a_l_poiss, b_l_poiss, a_alpha_bb, b_alpha_bb, a_s_bb, b_s_bb,
                                      tau_poiss,
                                      S_poiss, n_burnin_poiss, thin_poiss, seed)
  
  n_saved_iter_poiss <- length(output_poiss$lambda_vec)
  lambda_chain_poiss <- output_poiss$lambda_vec
  s_chain_poiss <- output_poiss$s_vec
  alpha_bar_chain_poiss <- output_poiss$alpha_bar_vec
  alpha_chain_poiss <- - alpha_bar_chain_poiss
  theta_chain_poiss <- s_chain_poiss + alpha_bar_chain_poiss
  
  ####### 4) Run BB + Negative-Binomial ################
  
  # Set tau for MALA
  tau_nb <- 0.001
  
  output_negbin <- gibbs_sampler_negbin_geometric(Z = train_mat,
                                                  nstar_0_nb, p_0_nb,  s_0_bb, alpha_bar_0_bb,
                                                  q_star_nb, alpha_p_nb, beta_p_nb, 
                                                  a_alpha_bb, b_alpha_bb, a_s_bb, b_s_bb,
                                                  tau_nb, fixed = c(F,F,F,F),
                                                  S_negbin, n_burnin_negbin, thin_negbin, seed)
  
  n_saved_iter_negbin <- length(output_negbin$nstar_vec)
  nstar_chain_negbin <- output_negbin$nstar_vec
  p_chain_negbin <- output_negbin$p_vec
  s_chain_negbin <- output_negbin$s_vec
  alpha_bar_chain_negbin <- output_negbin$alpha_bar_vec
  alpha_chain_negbin <- - alpha_bar_chain_negbin
  theta_chain_negbin <- s_chain_negbin + alpha_bar_chain_negbin
  
  ###### estimate limit distribution ################
  Kn = ncol(train_mat[,colSums(train_mat) > 0])
  
  ntilde_chain_poiss <- generate_Ntilde_chain_poiss(lambda_chain_poiss, alpha_chain_poiss, 
                                                    theta_chain_poiss, n = N, Kn)
  ntilde_chain_poiss <- as.data.frame(ntilde_chain_poiss)
  
  ntilde_chain_negbin <- generate_Ntilde_chain_negbin(nstar_chain_negbin, p_chain_negbin,
                                                      alpha_chain_negbin,
                                                      theta_chain_negbin, n = N, Kn)
  ntilde_chain_negbin <- as.data.frame(ntilde_chain_negbin)
  
  gg_list_poiss[[j]] <- ggplot(ntilde_chain_poiss, aes(x=ntilde_chain_poiss)) + 
    geom_density() + ggtitle(paste0("Poisson, xi = ", xi))
  gg_list_negbin[[j]] <- ggplot(ntilde_chain_negbin, aes(x=ntilde_chain_negbin)) + 
    geom_density() + ggtitle(paste0("Negbin, xi = ", xi))
  
  
  
}

ggarrange(gg_list_poiss[[1]], gg_list_poiss[[2]], gg_list_poiss[[3]],
          gg_list_poiss[[4]], gg_list_poiss[[5]], nrow = 2, ncol = 3)

ggarrange(gg_list_negbin[[1]], gg_list_negbin[[2]], gg_list_negbin[[3]],
          gg_list_negbin[[4]], gg_list_negbin[[5]], nrow = 2, ncol = 3)
