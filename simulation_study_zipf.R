####################################################################
#### Simulation using Zipf data - plotting the number of features -
## - EFPF and Method of moments & Poisson and Negative-Binomial
####################################################################

#################################################
######## Data Generation ########################
#################################################

seed = 1234

# total number of samples
#L = 30000
L = 8000

# value of xi
xi = 0.3

# generate data from zipf, in form of binary matrix
data_mat <- rzipf(n = L, K = 2*10^3, xi = xi, seed = seed)

# convert the binary matrix into list of features
data_list <- create_features_list(data_mat)
plot_trajectory(data_list)

print(paste0("Number of observed features: ", length(unique(unlist(data_list))) ))

# set the dimension of the training set and the training set itself
#Ns <- c(1000, 2000, 4000)
Ns <- c(1000, 2000, 4000)

############################################################
#### Fit the 4 models, under different training sets #######
############################################################

# list of ggplots - poiss
gg_list_poiss <- vector("list", length = length(Ns) )
# list of ggplots - negbin
gg_list_negbin <- vector("list", length = length(Ns) )

for (j in 1:length(Ns)){
  N <- Ns[j]
  
  train_list <- data_list[1:N]
  train_mat <- data_mat[1:N, ]
  
  # set the dimension of the test set and the test set itself
  M <- L - N
  test_list <- data_list[(N+1):L]
  test_mat <- data_mat[(N+1):L, ]
  train_counts <- tabulate(unlist(train_list))[tabulate(unlist(train_list)) != 0]
  
  ################################################################
  ################# Poisson model ################################
  ################################################################
  
  # EB: fix the model hyperparameters based on EFPF-maximization
  eb_efpf_poiss_BB <- EB_EFPF_poiss_BB(n = N, 
                                       counts = train_counts, pars_0 = c(-2, 5, 10))
  alpha_est_poiss_efpf <- eb_efpf_poiss_BB[1]
  theta_est_poiss_efpf <- eb_efpf_poiss_BB[2]
  lambda_est_poiss_efpf <- eb_efpf_poiss_BB[3]
  
  # EB: fix the model hyperparameters based on MM
  eb_mm_poiss_BB <- EB_MM_poiss_BB(n = N, 
                                   ntrain = round(N*2/3),
                                   val_rep = 5, data_list = train_list,
                                   pars_0 = c(-1,10,10), seed_id = 1234)
  alpha_est_poiss_mm <- eb_mm_poiss_BB[1]
  theta_est_poiss_mm <- eb_mm_poiss_BB[2]
  lambda_est_poiss_mm <- eb_mm_poiss_BB[3]
  
  # Compute E[KMN|Z_N] together with 90% credible intervals, for given N and train set
  CI_given_sample_poiss_efpf <- CI_Kmn_poiss_BB(alpha_est_poiss_efpf, theta_est_poiss_efpf,
                                          M, N, lambda_est_poiss_efpf, 0.95)
  CI_given_sample_poiss_mm <- CI_Kmn_poiss_BB(alpha_est_poiss_mm, theta_est_poiss_mm, M, N,
                                        lambda_est_poiss_mm, 0.95)
  gg_list_poiss[[j]] <- plot_Kmn_given_sample_with_observed_both(N, data_list, CI_given_sample_poiss_efpf, CI_given_sample_poiss_mm)
  
  
  
  ################################################################################
  ###########à Negative Binomial model ##########################################
  ##############################################################################
  
  # compute number of observed feature in the training set Kn
  Kn <- length(unique(unlist(train_list)))
  
  # EB: fix the model hyperparameters based on EFPF-maximization
  eb_efpf_negbin_BB <- EB_EFPF_fixed_negbin_mv_BB(n = N, 
                                                  counts = train_counts, pars_0 = c(-2, 5, 100, 150))
  
  alpha_est_negbin_efpf <- eb_efpf_negbin_BB[1]
  theta_est_negbin_efpf <- eb_efpf_negbin_BB[2]
  mu_est_negbin_efpf <- eb_efpf_negbin_BB[3]
  sigma2_est_negbin_efpf <- eb_efpf_negbin_BB[4]
  
  nstar_est_negbin_efpf <- (mu_est_negbin_efpf**2)/(sigma2_est_negbin_efpf - mu_est_negbin_efpf)
  p_est_negbin_efpf <- mu_est_negbin_efpf / sigma2_est_negbin_efpf
  

  # EB: fix the model hyperparameters based on MM
  eb_mm_negbin_mv_BB <- EB_MM_negbin_mv_BB(n = N, 
                                           ntrain = round(N*2/3),
                                           val_rep = 5, data_list = data_list,
                                           pars_0 = c(-2,100,100,5000), seed_id = 123)
  
  alpha_est_negbin_mm <- eb_mm_negbin_mv_BB[1]
  theta_est_negbin_mm <- eb_mm_negbin_mv_BB[2]
  mu_est_negbin_mm <- eb_mm_negbin_mv_BB[3]
  sigma2_est_negbin_mm <- eb_mm_negbin_mv_BB[4]
  
  nstar_est_negbin_mm <- (mu_est_negbin_mm**2)/(sigma2_est_negbin_mm - mu_est_negbin_mm)
  p_est_negbin_mm <- mu_est_negbin_mm / sigma2_est_negbin_mm
  
  # Compute E[KMN|Z_N] together with 90% credible intervals, for given N and train set
  CI_given_sample_negbin_efpf <- CI_Kmn_negbin_BB(alpha_est_negbin_efpf, theta_est_negbin_efpf,
                                                  M, N, Kn, nstar_est_negbin_efpf, 
                                                  p_est_negbin_efpf, 0.95)
  CI_given_sample_negbin_mm <- CI_Kmn_negbin_BB(alpha_est_negbin_mm, theta_est_negbin_mm,
                                                  M, N, Kn, nstar_est_negbin_mm, 
                                                  p_est_negbin_mm, 0.95)
  gg_list_negbin[[j]] <- plot_Kmn_given_sample_with_observed_both(N, data_list, CI_given_sample_negbin_efpf, CI_given_sample_negbin_mm)
  
  
}

#save(gg_list_poiss, file = "zipf_gg_poiss")
#save(gg_list_negbin, file = "zipf_gg_negbin")
load("zipf_gg_poiss")
load("zipf_gg_negbin")

# display plots on the same page
ggarrange(gg_list_poiss[[1]], gg_list_poiss[[2]], gg_list_poiss[[3]],
          gg_list_negbin[[1]], gg_list_negbin[[2]], gg_list_negbin[[3]],
          ncol = 3, nrow = 2)



####################################################################
#### Simulation using Zipf data - plotting the number of features -
## - EFPF with Poisson and Negative-Binomial - same setting
####################################################################

#################################################
######## Data Generation ########################
#################################################

seed = 1234

# total number of samples
#L = 30000
L = 8000

# value of xi
xi = 0.3

# generate data from zipf, in form of binary matrix
data_mat <- rzipf(n = L, K = 2*10^3, xi = xi, seed = seed)

# convert the binary matrix into list of features
data_list <- create_features_list(data_mat)
plot_trajectory(data_list)

print(paste0("Number of observed features: ", length(unique(unlist(data_list))) ))

# set the dimension of the training set and the training set itself
#Ns <- c(1000, 2000, 4000)
Ns <- c(1000, 2000, 4000)

############################################################
#### Fit the 4 models, under different training sets #######
############################################################

# list of ggplots - poiss
gg_list_poiss_efpf <- vector("list", length = length(Ns) )
# list of ggplots - negbin
gg_list_negbin_efpf <- vector("list", length = length(Ns) )

for (j in 1:length(Ns)){
  N <- Ns[j]
  
  train_list <- data_list[1:N]
  train_mat <- data_mat[1:N, ]
  
  # set the dimension of the test set and the test set itself
  M <- L - N
  test_list <- data_list[(N+1):L]
  test_mat <- data_mat[(N+1):L, ]
  train_counts <- tabulate(unlist(train_list))[tabulate(unlist(train_list)) != 0]
  
  ################################################################
  ################# Poisson model ################################
  ################################################################
  
  # EB: fix the model hyperparameters based on EFPF-maximization
  eb_efpf_poiss_BB <- EB_EFPF_poiss_BB(n = N, 
                                       counts = train_counts, pars_0 = c(-2, 5, 10))
  alpha_est_poiss_efpf <- eb_efpf_poiss_BB[1]
  theta_est_poiss_efpf <- eb_efpf_poiss_BB[2]
  lambda_est_poiss_efpf <- eb_efpf_poiss_BB[3]
  
  # Compute E[KMN|Z_N] together with 90% credible intervals, for given N and train set
  CI_given_sample_poiss_efpf <- CI_Kmn_poiss_BB(alpha_est_poiss_efpf, theta_est_poiss_efpf,
                                                M, N, lambda_est_poiss_efpf, 0.95)
  gg_list_poiss_efpf[[j]] <- plot_Kmn_given_sample_with_observed(N, data_list, CI_given_sample_poiss_efpf)
  
  
  ################################################################################
  ###########à Negative Binomial model ##########################################
  ##############################################################################
  
  # compute number of observed feature in the training set Kn
  Kn <- length(unique(unlist(train_list)))
  
  # EB: fix the model hyperparameters based on EFPF-maximization
  eb_efpf_negbin_BB <- EB_EFPF_fixed_negbin_mv_BB(n = N, 
                                                  counts = train_counts, pars_0 = c(-2, 5, 100, 150))
  
  alpha_est_negbin_efpf <- eb_efpf_negbin_BB[1]
  theta_est_negbin_efpf <- eb_efpf_negbin_BB[2]
  mu_est_negbin_efpf <- eb_efpf_negbin_BB[3]
  sigma2_est_negbin_efpf <- eb_efpf_negbin_BB[4]
  
  nstar_est_negbin_efpf <- (mu_est_negbin_efpf**2)/(sigma2_est_negbin_efpf - mu_est_negbin_efpf)
  p_est_negbin_efpf <- mu_est_negbin_efpf / sigma2_est_negbin_efpf
  
  
  # Compute E[KMN|Z_N] together with 90% credible intervals, for given N and train set
  CI_given_sample_negbin_efpf <- CI_Kmn_negbin_BB(alpha_est_negbin_efpf, theta_est_negbin_efpf,
                                                  M, N, Kn, nstar_est_negbin_efpf, 
                                                  p_est_negbin_efpf, 0.95)
  gg_list_negbin_efpf[[j]] <- plot_Kmn_given_sample_with_observed(N, data_list, CI_given_sample_negbin_efpf)
  
  
}

#save(gg_list_poiss_efpf, file = "zipf_gg_poiss_efpf")
#save(gg_list_negbin_efpf, file = "zipf_gg_negbin_efpf")
load("zipf_gg_poiss_efpf")
load("zipf_gg_negbin_efpf")

# display plots on the same page
ggarrange(gg_list_poiss_efpf[[1]], gg_list_poiss_efpf[[2]], gg_list_poiss_efpf[[3]],
          gg_list_negbin_efpf[[1]], gg_list_negbin_efpf[[2]], gg_list_negbin_efpf[[3]],
          ncol = 3, nrow = 2)



####################################################################
## Simulation using Zipf data - accuracy over replicated datasets ##
######### Compare Poisson and Negative-Binomial models #############
#################### using EFPF-maximization #######################
####################################################################
library(scales)

seed = 1234

# total number of samples
L = 8000

# value of xi
xis = c(0.3, 0.6, 0.8, 1)

# total number of datasets
S = 100

# set the dimension of the training set and the training set itself
Ns <- c(50, 200, 1000)

boxplot_acc <- data.frame(matrix(ncol = 4, nrow = 0))
coln <- c("N", "acc", "Model", "Xi")
colnames(boxplot_acc) <- coln

for (h in 1:length(xis)){
  xi <- xis[h]
  # matrix of accuracy values, dimensions: length(Ns) x S 
  perc_acc_mat_poiss <- matrix(NA, nrow = length(Ns), ncol = S)
  perc_acc_mat_negbin <- matrix(NA, nrow = length(Ns), ncol = S)
  
  s = 1
  while (s <= S){
    s <- s + 1
    # generate data from zipf, in form of binary matrix
    data_mat <- rzipf(n = L, K = 2*10^3, xi = xi, seed = seed)
    
    # convert the binary matrix into list of features
    data_list <- create_features_list(data_mat)
    
    for (j in 1:length(Ns)){
      N <- Ns[j]
      
      train_list <- data_list[1:N]
      train_mat <- data_mat[1:N, ]
      
      # set the dimension of the test set and the test set itself
      M <- L - N
      test_list <- data_list[(N+1):L]
      test_mat <- data_mat[(N+1):L, ]
      
      # compute number of observed feature in the training set Kn
      Kn <- length(unique(unlist(train_list)))
      
      # EB: fix the model hyperparameters based on EFPF-maximization
      train_counts <- colSums(train_mat)[colSums(train_mat)!=0]
      
      # Fit Poisson model
      eb_efpf_poiss_BB <- EB_EFPF_poiss_BB(n = N, 
                                           counts = train_counts, pars_0 = c(-5, 10, 100))
      
      alpha_est_poiss <- eb_efpf_poiss_BB[1]
      theta_est_poiss <- eb_efpf_poiss_BB[2]
      lambda_est <- eb_efpf_poiss_BB[3]
      
      # Fit Negative-Binomial model
      eb_efpf_negbin_BB <- EB_EFPF_fixed_negbin_mv_BB(n = N, 
                                                      counts = train_counts, pars_0 = c(-2, 5, 100, 150))
      
      alpha_est_negbin <- eb_efpf_negbin_BB[1]
      theta_est_negbin <- eb_efpf_negbin_BB[2]
      mu_est <- eb_efpf_negbin_BB[3]
      sigma2_est <- eb_efpf_negbin_BB[4]
      
      nstar_est <- (mu_est**2)/(sigma2_est - mu_est)
      p_est <- mu_est / sigma2_est
      
      ############################### FOR ACCURACY  ##################################
      # estimated number of new features for poisson
      est_new_features_poiss <- mean_kmn_poiss_BB(alpha = alpha_est_poiss, theta = theta_est_poiss, 
                                                  m = M, n = N, lambda = lambda_est)
      # estimated number of new features for negative-binomial
      est_new_features_negbin <- mean_kmn_negbin_BB(alpha = alpha_est_negbin, theta = theta_est_negbin, 
                                                    m = M, n = N, Kn = Kn, nstar = nstar_est, p = p_est)
      
      if (is.nan(est_new_features_poiss) | is.nan(est_new_features_negbin)){
        s <- s - 1
        break
      }
      else{
        # compute the percentage accuracy of the estimate with respect to the test set - poisson
        perc_acc_mat_poiss[j,s-1] <- perc_accuracy(train_list, test_list, est_new_features_poiss)
        # compute the percentage accuracy of the estimate with respect to the test set - negative binomial
        perc_acc_mat_negbin[j,s-1] <- perc_accuracy(train_list, test_list, est_new_features_negbin)
      }  
    }
    
    seed = seed + 1
    
  }
  
  box_data_poiss <- data.frame(
    N = as.factor(rep(Ns, each = S)),
    acc = c(t(perc_acc_mat_poiss)),
    Model = rep("Poiss-BB", S*length(Ns))
  )
  
  box_data_negbin <- data.frame(
    N = as.factor(rep(Ns, each = S)),
    acc = c(t(perc_acc_mat_negbin)),
    Model = rep("NB-BB", S*length(Ns))
  )
  
  box_data <- rbind(box_data_poiss, box_data_negbin)
  box_data$Model <- factor(box_data$Model, levels = c("Poiss-BB", "NB-BB"))
  box_data$Xi <- rep(xi, 2*length(Ns)*S)
  
  boxplot_acc <- rbind(boxplot_acc, box_data)
  
}


# plots
ggplot(boxplot_acc, aes(x=N, y=acc, fill=Model)) + 
  geom_boxplot() + 
  facet_wrap(~Xi) +
  xlab("# training set") + ylab("Accuracy") + theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(paste0("Zipf: average on S = ", S)) +
  scale_y_continuous(breaks = pretty_breaks()) 

