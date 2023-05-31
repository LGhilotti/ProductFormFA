####
##### MODEL 2: the random uniform model ############################
####

rm(list=ls())
library(ggmcmc)
library(coda)
library(ggpubr)
library(scales)

##### Single dataset -> Ntilde (Poiss/NB) and extrapolation (Poiss/NB/Gamma) #####

###### 1) Read results:  MCMC convergence ####################
load(file = "chao_model_simulation/m2/m2_params_poiss.Rda")
load(file =  "chao_model_simulation/m2/m2_params_negbin.Rda")
load(file =  "chao_model_simulation/m2/m2_params_ibp.Rda", )

###### 2) Read results: samples from limiting distributions (Poiss/NB) ##############
load(file = "chao_model_simulation/m2/m2_ntilde_poiss.Rda")
load(file = "chao_model_simulation/m2/m2_ntilde_negbin.Rda")

###### 3) Read results: CI for extrapolation (Poiss/NB/Gamma) ################
list_kmn_pred_test_poiss <- readRDS(file = "chao_model_simulation/m2/m2_ci_poiss.rds")
list_kmn_pred_test_negbin <- readRDS(file = "chao_model_simulation/m2/m2_ci_negbin.rds")
list_kmn_pred_test_ibp <- readRDS(file = "chao_model_simulation/m2/m2_ci_ibp.rds")

###### 4) Read the data ###############################
data_mat <- readRDS(file = "chao_model_simulation/m2/m2_data_mat.rds")
L <- nrow(data_mat)
Ms <- sapply(list_kmn_pred_test_poiss, function(l) length(l$medians))
Ns <- L - Ms


###### 5) Check MCMC convergence###################

###### check mcmc mixing poisson
samples_poiss <- mcmc.list(mcmc(params_poiss))
samples_ggs_poiss <- ggs(samples_poiss, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_poiss) + 
  facet_wrap(~Parameter, nrow = length(Ns), scales = "free")

###### check mcmc mixing negbin
samples_negbin <- mcmc.list(mcmc(params_negbin))
samples_ggs_negbin <- ggs(samples_negbin, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_negbin) + 
  facet_wrap(~Parameter, nrow = length(Ns), scales = "free")

###### check mcmc mixing Gamma ibp
samples_ibp <- mcmc.list(mcmc(params_ibp))
samples_ggs_ibp <- ggs(samples_ibp, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_ibp) + 
  facet_wrap(~Parameter, nrow = length(Ns), scales = "free")


####### 6) Plot limiting distributions (Poiss/NB) ##############
# Poisson
samples_ntilde_poiss <- mcmc.list(mcmc(gg_ntilde_poiss))
samples_ntilde_ggs_poiss <- ggs(samples_ntilde_poiss, keep_original_order = TRUE)
ggs_density(samples_ntilde_ggs_poiss) + 
  facet_wrap(~Parameter, nrow = length(Ns), scales = "free")

# Negative Binomial
samples_ntilde_negbin <- mcmc.list(mcmc(gg_ntilde_negbin))
samples_ntilde_ggs_negbin <- ggs(samples_ntilde_negbin, keep_original_order = TRUE)
ggs_density(samples_ntilde_ggs_negbin) + 
  facet_wrap(~Parameter, nrow = length(Ns), scales = "free")

######## 7) Plot Extrapolation curve (Poiss/NB/Gamma) ################
gg_kmn_pred_test_poiss  <- vector(mode="list", length = length(Ns))
gg_kmn_pred_test_negbin  <- vector(mode="list", length = length(Ns))
gg_kmn_pred_test_ibp  <- vector(mode="list", length = length(Ns))

for (j in 1:length(Ns)){
  N <- Ns[j]
  M <- L - N
  
  train_mat <- data_mat[1:N,]
  test_mat <- data_mat[(N+1):L, ]
  # convert the binary matrix into list of features
  train_list <- create_features_list(train_mat)
  test_list <- create_features_list(test_mat)
  
  
  gg_kmn_pred_test_poiss_ <- plot_Kmn_median_pred_and_test(train_list = train_list, test_list = test_list,
                                                           ci = list_kmn_pred_test_poiss[[paste0("N.",N)]], n_avg = 100)
  gg_kmn_pred_test_poiss[[j]] <- gg_kmn_pred_test_poiss_ + ggtitle(paste0("Poiss, N = ", N))
  
  gg_kmn_pred_test_negbin_ <- plot_Kmn_median_pred_and_test(train_list = train_list, test_list = test_list,
                                                            ci = list_kmn_pred_test_negbin[[paste0("N.",N)]], n_avg = 100)
  gg_kmn_pred_test_negbin[[j]] <- gg_kmn_pred_test_negbin_ + ggtitle(paste0("NegBin, N = ", N))
  
  gg_kmn_pred_test_ibp_ <- plot_Kmn_median_pred_and_test(train_list = train_list, test_list = test_list,
                                                         ci = list_kmn_pred_test_ibp[[paste0("N.",N)]], n_avg = 100)
  gg_kmn_pred_test_ibp[[j]] <- gg_kmn_pred_test_ibp_ + ggtitle(paste0("Gamma IBP, N = ", N))
  
}

# Poisson
ggarrange(plotlist = gg_kmn_pred_test_poiss)
ggsave("chao_model_simulation/m2/m2_poiss_extr.pdf", width = 16, height = 8, units = "in")

# Negative Binomial
ggarrange(plotlist = gg_kmn_pred_test_negbin)
ggsave("chao_model_simulation/m2/m2_negbin_extr.pdf", width = 16, height = 8, units = "in")

# Gamma IBP
ggarrange(plotlist = gg_kmn_pred_test_ibp)
ggsave("chao_model_simulation/m2/m2_ibp_extr.pdf", width = 16, height = 8, units = "in")





##### Accuracy on multiple datasets #####

###### 1) Read results: limit distribution estimates #####
load(file = "chao_model_simulation/m2/m2_avg_ntilde_poiss.Rda")
load(file = "chao_model_simulation/m2/m2_avg_ntilde_negbin.Rda")

###### 2) Read results: boxplots on accuracy #####
load(file = "chao_model_simulation/m2/m2_acc_poiss.Rda")
load(file = "chao_model_simulation/m2/m2_acc_negbin.Rda")
load(file = "chao_model_simulation/m2/m2_acc_ibp.Rda")

D <- nrow(avg_ntilde_poiss)

###### 3) Plot limit distribution estimate #####
avg_ntilde_poiss_long <- gather(avg_ntilde_poiss, training, estimate, 
                                paste0("N.",Ns[1]):paste0("N.",Ns[length(Ns)]), factor_key=TRUE) %>%
  add_column(model = "Poiss", N = rep(Ns[1:length(Ns)], each = D))

avg_ntilde_negbin_long <- gather(avg_ntilde_negbin, training, estimate, 
                                 paste0("N.",Ns[1]):paste0("N.",Ns[length(Ns)]), factor_key=TRUE) %>%
  add_column(model = "Neg-Bin", N = rep(Ns[1:length(Ns)], each = D))

joint_ntilde_long <- bind_rows(avg_ntilde_poiss_long, avg_ntilde_negbin_long)

# plots
ggplot(joint_ntilde_long, aes( y=estimate, fill=model)) + 
  geom_boxplot() + 
  facet_wrap(~N) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(paste0("Model 1: average on D = ", D)) +
  scale_y_continuous(breaks = pretty_breaks()) 



###### 4) Plot boxplots on accuracy #####
acc_poiss_long <- gather(acc_poiss, training, accuracy, 
                         paste0("N.",Ns[1]):paste0("N.",Ns[length(Ns)]), factor_key=TRUE) %>%
  add_column(model = "Poiss", N = rep(Ns, each = D))

acc_negbin_long <- gather(acc_negbin, training , accuracy, 
                          paste0("N.",Ns[1]):paste0("N.",Ns[length(Ns)]), factor_key=TRUE) %>%
  add_column(model = "Neg-Bin", N = rep(Ns, each = D))

acc_ibp_long <- gather(acc_ibp, training , accuracy, 
                       paste0("N.",Ns[1]):paste0("N.",Ns[length(Ns)]), factor_key=TRUE) %>%
  add_column(model = "Gamma IBP", N = rep(Ns, each = D))

joint_long <- bind_rows(acc_poiss_long, acc_negbin_long, acc_ibp_long)

# plots
ggplot(joint_long, aes( y=accuracy, fill=model)) + 
  geom_boxplot() + 
  facet_wrap(~N) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(paste0("Model 1: average on D = ", D)) +
  scale_y_continuous(breaks = pretty_breaks()) 

