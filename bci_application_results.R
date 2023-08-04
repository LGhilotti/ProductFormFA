##############################################################
############ REAL DATA APPLICATION: BCI (vegan) ###################
#############################################################
rm(list=ls())

library(ProductFormFA)

library(ggmcmc)
library(coda)
library(ggpubr)
library(scales)
library(tidyverse)
library(gridExtra)
library(grid)
library(patchwork)

##########################################################
#### (A) Training/test approach ############################
##########################################################

###### 1) Read results:  MCMC convergence ####################
load(file = "bci_results/bci_params_poiss.Rda")
load(file =  "bci_results/bci_params_negbin.Rda")
load(file =  "bci_results/bci_params_ibp.Rda" )
load(file =  "bci_results/bci_params_sp.Rda" )

###### 2) Read results: samples from limiting distributions (Poiss/NB) ##############
load(file = "bci_results/bci_ntilde_poiss.Rda")
load(file = "bci_results/bci_ntilde_negbin.Rda")

###### 3) Read results: CI for extrapolation (Poiss/NB/Gamma) ################
list_kmn_pred_test_poiss <- readRDS(file = "bci_results/bci_ci_poiss.rds")
list_kmn_pred_test_negbin <- readRDS(file = "bci_results/bci_ci_negbin.rds")
list_kmn_pred_test_ibp <- readRDS(file = "bci_results/bci_ci_ibp.rds")
list_kmn_pred_test_sp <- readRDS(file = "bci_results/bci_ci_sp.rds")

###### 3.b) Read results: CI for insample (Poiss/NB/Gamma) ################
list_kn_rarefaction_poiss <- readRDS(file = "bci_results/bci_ci_insample_poiss.rds")
list_kn_rarefaction_negbin <- readRDS(file = "bci_results/bci_ci_insample_negbin.rds")
list_kn_rarefaction_ibp <- readRDS(file = "bci_results/bci_ci_insample_ibp.rds")
list_kn_rarefaction_sp <- readRDS(file = "bci_results/bci_ci_insample_sp.rds")


###### 4) Read the data ###############################
data_mat <- readRDS(file = "bci_results/bci_data_mat.rds")
data_list <- create_features_list(data_mat)
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

###### check mcmc mixing SB-SP
samples_sp <- mcmc.list(mcmc(params_sp))
samples_ggs_sp <- ggs(samples_sp, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_sp) + 
  facet_wrap(~Parameter, nrow = length(Ns), scales = "free")


###### 6) Rarefaction curves
gg_kn_rarefaction_all  <- vector(mode="list", length = length(Ns))

for (j in 1:length(Ns)){
  N <- Ns[j]
  
  train_mat <- data_mat[1:N,]
  # convert the binary matrix into list of features
  train_list <- create_features_list(train_mat)
  
  gg_kn_rarefaction_all[[j]] <- plot_Kn_median_and_rarefaction_all(
    train_list = train_list,
    ci_poiss = list_kn_rarefaction_poiss[[paste0("N.",N)]], 
    ci_negbin = list_kn_rarefaction_negbin[[paste0("N.",N)]],
    ci_ibp = list_kn_rarefaction_ibp[[paste0("N.",N)]],
    n_avg = 100) +
    ggtitle(paste0("n = ", N)) + theme(plot.title = element_text(size=12)) 
  
  
  if (j != 1){
    gg_kn_rarefaction_all[[j]] <- gg_kn_rarefaction_all[[j]] +
      theme(#axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        axis.title.y = element_blank() )
  }
  
}

# Print plots
fig_bci_rare_ <- wrap_plots(gg_kn_rarefaction_all, nrow = 1, ncol = 3) +   plot_layout(guides = "collect") & theme(legend.position = 'right') & xlab(NULL) & theme(plot.margin = margin(5.5, 5.5, 0, 5.5))

fig_bci_rare <- wrap_elements(panel = fig_bci_rare_) +
  labs(tag = "# observations") +
  theme(
    plot.tag = element_text(size = rel(1)),
    plot.tag.position = "bottom") 

ggsave(filename = "Plots_paper/plot_bci_rare.png", width = 10, height = 4, dpi = 300, units = "in", device='png')


###### 7) Prediction
gg_kmn_pred_test_all  <- vector(mode="list", length = length(Ns))

for (j in 1:length(Ns)){
  N <- Ns[j]
  M <- L - N
  
  train_mat <- data_mat[1:N,]
  test_mat <- data_mat[(N+1):L, ]
  # convert the binary matrix into list of features
  train_list <- create_features_list(train_mat)
  test_list <- create_features_list(test_mat)
  
  
  gg_kmn_pred_test_all[[j]] <- plot_Kmn_median_pred_and_test_all(
    train_list = train_list,
    test_list = test_list,
    ci_poiss = list_kmn_pred_test_poiss[[paste0("N.",N)]], 
    ci_negbin = list_kmn_pred_test_negbin[[paste0("N.",N)]],
    ci_ibp = list_kmn_pred_test_ibp[[paste0("N.",N)]],
    n_avg = 100) +
    ggtitle(paste0("n = ", N)) + theme(plot.title = element_text(size=12))
  
  
  if (j != 1){
    gg_kmn_pred_test_all[[j]] <- gg_kmn_pred_test_all[[j]] +
      theme(#axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        axis.title.y = element_blank() )
  }
  
}

# Print plots
fig_bci_pred_ <- wrap_plots(gg_kmn_pred_test_all, nrow = 1, ncol = 3) +   plot_layout(guides = "collect") & theme(legend.position = 'right') & xlab(NULL) & theme(plot.margin = margin(5.5, 5.5, 0, 5.5))

fig_bci_pred <- wrap_elements(panel = fig_bci_pred_) +
  labs(tag = "# observations") +
  theme(
    plot.tag = element_text(size = rel(1)),
    plot.tag.position = "bottom") 

ggsave(filename = "Plots_paper/plot_bci_pred.png", width = 10, height = 4, dpi = 300, units = "in", device='png')


###### 8) Richness
gg_ntilde_poiss_long <- gather(gg_ntilde_poiss, training, estimate, 
                               paste0("N.",Ns[1]):paste0("N.",Ns[length(Ns)]), 
                               factor_key=TRUE) %>%
  add_column(Model= "BBmixP", N = rep(Ns, each = nrow(gg_ntilde_poiss)))


gg_ntilde_negbin_long <- gather(gg_ntilde_negbin, training, estimate, 
                                paste0("N.",Ns[1]):paste0("N.",Ns[length(Ns)]), 
                                factor_key=TRUE) %>%
  add_column(Model= "BBmixNB", N = rep(Ns, each = nrow(gg_ntilde_negbin)))


joint_total_long <- bind_rows(gg_ntilde_poiss_long, gg_ntilde_negbin_long) %>%
  mutate(Model = fct_relevel(Model, c("BBmixP", "BBmixNB"))) 

# plot
ggplot(joint_total_long, aes(x = estimate, color = Model)) + 
  #geom_density(alpha = 0.5, linewidth = 0.8) +
  stat_density(aes(x=estimate, colour=Model),
               geom="line",position="identity") +
  geom_vline(xintercept = ncol(data_mat), color="grey", linetype="dashed", linewidth=0.8) + # observed features in the whole dataset
  facet_wrap(~N, labeller = labeller(N = ~ paste("n = ", .x)), scales = "free", nrow = 1) +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 12))+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(breaks = pretty_breaks()) +
  xlab("# distinct features") + rremove("ylab") +
  scale_color_manual(values = c("BBmixP" = "forestgreen",
                                "BBmixNB" = "royalblue1")) 

ggsave(filename = "Plots_paper/plot_bci_richness.png", width = 10, height = 4, dpi = 300, units = "in", device='png')



##########################################################
#### (B) Only prediction #################################
##########################################################
rm(list=ls())
###### 1) Read results:  MCMC convergence ####################
load(file = "bci_results/bci_op_params_poiss.Rda")
load(file =  "bci_results/bci_op_params_negbin.Rda")
load(file =  "bci_results/bci_op_params_ibp.Rda" )
load(file =  "bci_results/bci_op_params_sp.Rda" )

lambda_chain_poiss <- params_poiss[["lambda"]]  
alpha_chain_poiss <- params_poiss[["alpha"]]
theta_chain_poiss <- params_poiss[["theta"]]

nstar_chain_negbin <- params_negbin[["nstar"]] 
p_chain_negbin <- params_negbin[["p"]]
alpha_chain_negbin <- params_negbin[["alpha"]]
theta_chain_negbin <- params_negbin[["theta"]]

a_chain_ibp <- params_ibp[["a"]]
b_chain_ibp <- params_ibp[["b"]]
alpha_chain_ibp <- params_ibp[["alpha"]]
theta_chain_ibp <- params_ibp[["theta"]] 

c_chain_sp <- params_sp[["c"]]
beta_chain_sp <- params_sp[["beta"]] 
alpha_chain_sp <- params_sp[["alpha"]]

###### 2) Read results: samples from limiting distributions (Poiss/NB) ##############
load(file = "bci_results/bci_op_ntilde_poiss.Rda")
load(file = "bci_results/bci_op_ntilde_negbin.Rda")

###### 3.b) Read results: CI for insample (Poiss/NB/Gamma) ################
kn_rarefaction_poiss <- readRDS(file = "bci_results/bci_op_ci_insample_poiss.rds")
kn_rarefaction_negbin <- readRDS(file = "bci_results/bci_op_ci_insample_negbin.rds")
kn_rarefaction_ibp <- readRDS(file = "bci_results/bci_op_ci_insample_ibp.rds")
kn_rarefaction_sp <- readRDS(file = "bci_results/bci_op_ci_insample_sp.rds")

###### 4) Read the data ###############################
data_mat <- readRDS(file = "bci_results/bci_data_mat.rds")
data_list <- create_features_list(data_mat)
L <- nrow(data_mat)


###### 5) Check MCMC convergence###################

###### check mcmc mixing poisson
samples_poiss <- mcmc.list(mcmc(params_poiss))
samples_ggs_poiss <- ggs(samples_poiss, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_poiss) 

###### check mcmc mixing negbin
samples_negbin <- mcmc.list(mcmc(params_negbin))
samples_ggs_negbin <- ggs(samples_negbin, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_negbin) 

###### check mcmc mixing Gamma ibp
samples_ibp <- mcmc.list(mcmc(params_ibp))
samples_ggs_ibp <- ggs(samples_ibp, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_ibp) 

###### check mcmc mixing SB-SP
samples_sp <- mcmc.list(mcmc(params_sp))
samples_ggs_sp <- ggs(samples_sp, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_sp) 

###### 6) Rarefaction curves
gg_kn_rarefaction_all <- plot_Kn_median_and_rarefaction_all(
  train_list = data_list,
  ci_poiss = kn_rarefaction_poiss, 
  ci_negbin = kn_rarefaction_negbin,
  ci_ibp = kn_rarefaction_ibp,
  n_avg = 100) +
  ggtitle(paste0("n = ", L)) + theme(plot.title = element_text(size=12)) 

ggsave(filename = "Plots_paper/plot_bci_op_rare.png", width = 10, height = 4, dpi = 300, units = "in", device='png')


# ###### 7.a) Prediction on different horizons (computation) - 1000 new subjects ##########
# 
# Ts <- c(1000) # horizon ahead to check prediction
# 
# list_kmn_pred_poiss <- vector(mode="list", length = length(Ts))
# names(list_kmn_pred_poiss) <- paste("T", Ts, sep = ".")
# list_kmn_pred_negbin <- vector(mode="list", length = length(Ts))
# names(list_kmn_pred_negbin) <- paste("T", Ts, sep = ".")
# list_kmn_pred_ibp <- vector(mode="list", length = length(Ts))
# names(list_kmn_pred_ibp) <- paste("T", Ts, sep = ".")
# list_kmn_pred_sp <- vector(mode="list", length = length(Ts))
# names(list_kmn_pred_sp) <- paste("T", Ts, sep = ".")
# 
# Kn = ncol(data_mat[,colSums(data_mat) > 0])
# 
# for (j in 1:length(Ts)){
#   hor <- Ts[j]
#   
#   # Poisson
#   kmn_chain_poiss <- generate_Kmn_chain_poiss(lambda_chain_poiss, alpha_chain_poiss,
#                                               theta_chain_poiss, M = hor, n = L)
#   
#   est_ci_pred_poiss <- matrix(NA, nrow = hor, ncol = 3)
#   # first column = lower bound
#   # second columns = medians
#   # third columns = upper bound
#   for (m in 1:hor){
#     est_ci_pred_poiss[m,] <- quantile(kmn_chain_poiss[m,], probs = c(0.025,0.5,0.975))
#   }
#   est_ci_pred_poiss <- list("medians" = est_ci_pred_poiss[,2],
#                             "lbs" = est_ci_pred_poiss[,1],
#                             "ubs" = est_ci_pred_poiss[,3])
#   
#   list_kmn_pred_poiss[[paste0("T.",hor)]] <- est_ci_pred_poiss
#   
#   # Negative Binomial
#   kmn_chain_negbin <- generate_Kmn_chain_negbin(nstar_chain_negbin, p_chain_negbin,
#                                                 alpha_chain_negbin, theta_chain_negbin,
#                                                 M = hor, n = L, Kn)
#   
#   est_ci_pred_negbin <- matrix(NA, nrow = hor, ncol = 3)
#   # first column = lower bound
#   # second columns = medians
#   # third columns = upper bound
#   for (m in 1:hor){
#     est_ci_pred_negbin[m,] <- quantile(kmn_chain_negbin[m,], probs = c(0.025,0.5,0.975))
#   }
#   est_ci_pred_negbin <- list("medians" = est_ci_pred_negbin[,2],
#                              "lbs" = est_ci_pred_negbin[,1],
#                              "ubs" = est_ci_pred_negbin[,3])
#   
#   list_kmn_pred_negbin[[paste0("T.",hor)]] <- est_ci_pred_negbin
#   
#   
#   # IBP + Gamma
#   kmn_chain_ibp <- generate_Kmn_chain_gamma_ibp(a_chain_ibp, b_chain_ibp, alpha_chain_ibp,
#                                                 theta_chain_ibp, M = hor, n = L, Kn)
#   
#   est_ci_pred_ibp <- matrix(NA, nrow = hor, ncol = 3)
#   # first column = lower bound
#   # second columns = medians
#   # third columns = upper bound
#   for (m in 1:hor){
#     est_ci_pred_ibp[m,] <- quantile(kmn_chain_ibp[m,], probs = c(0.025,0.5,0.975))
#   }
#   est_ci_pred_ibp <- list("medians" = est_ci_pred_ibp[,2],
#                           "lbs" = est_ci_pred_ibp[,1],
#                           "ubs" = est_ci_pred_ibp[,3])
#   
#   list_kmn_pred_ibp[[paste0("T.",hor)]] <- est_ci_pred_ibp
#   
#   # SB-SP
#   kmn_chain_sp <- generate_Kmn_chain_gamma_ibp(a_chain = c_chain_sp + 1,
#                                                b_chain = beta_chain_sp*(1-alpha_chain_sp)/alpha_chain_sp, 
#                                                alpha_chain = alpha_chain_sp,
#                                                theta_chain = 1 - alpha_chain_sp,
#                                                M = hor, n = L, Kn)
#   
#   est_ci_pred_sp <- matrix(NA, nrow = hor, ncol = 3)
#   # first column = lower bound
#   # second columns = medians
#   # third columns = upper bound
#   for (m in 1:hor){
#     est_ci_pred_sp[m,] <- quantile(kmn_chain_sp[m,], probs = c(0.025,0.5,0.975))
#   }
#   est_ci_pred_sp <- list("medians" = est_ci_pred_sp[,2],
#                          "lbs" = est_ci_pred_sp[,1],
#                          "ubs" = est_ci_pred_sp[,3])
#   
#   list_kmn_pred_sp[[paste0("T.",hor)]] <- est_ci_pred_sp
#   
#   
# }
# 
# ######### Save results: CI for extrapolation (Poiss/NB/Gamma/SP)
# # # Poisson
# saveRDS(list_kmn_pred_poiss, "bci_results/bci_op_1000_ci_poiss.rds")
# # # Negative Binomial
# saveRDS(list_kmn_pred_negbin, "bci_results/bci_op_1000_ci_negbin.rds")
# # # Gamma IBP
# saveRDS(list_kmn_pred_ibp, "bci_results/bci_op_1000_ci_ibp.rds")
# # # SB-SP
# saveRDS(list_kmn_pred_sp, "bci_results/bci_op_1000_ci_sp.rds")


###### 7.b) Prediction on different horizons (results) - 1000 new subjects ########
list_kmn_pred_poiss <- readRDS(file = "bci_results/bci_op_1000_ci_poiss.rds")
list_kmn_pred_negbin <- readRDS(file = "bci_results/bci_op_1000_ci_negbin.rds")
list_kmn_pred_ibp <- readRDS(file = "bci_results/bci_op_1000_ci_ibp.rds")
list_kmn_pred_sp <- readRDS(file = "bci_results/bci_op_1000_ci_sp.rds")

hor <- 1000 # horizon ahead to check prediction

gg_kmn_pred_all <- plot_Kmn_median_pred_all(
  data_list = data_list,
  ci_poiss = list_kmn_pred_poiss[[paste0("T.",hor)]], 
  ci_negbin = list_kmn_pred_negbin[[paste0("T.",hor)]],
  ci_ibp = list_kmn_pred_ibp[[paste0("T.",hor)]],
  n_avg = 100) + ggtitle(paste0("# new observations = ", hor))

ggsave(filename = "Plots_paper/plot_bci_op_1000_pred.png", width = 8, height = 6, dpi = 300, units = "in", device='png')


####### 7.c) Prediction on different horizons (computation) - small-set of new subjects ######

# Ts <- c(200, 600) # horizon ahead to check prediction
# 
# list_kmn_pred_poiss <- vector(mode="list", length = length(Ts))
# names(list_kmn_pred_poiss) <- paste("T", Ts, sep = ".")
# list_kmn_pred_negbin <- vector(mode="list", length = length(Ts))
# names(list_kmn_pred_negbin) <- paste("T", Ts, sep = ".")
# list_kmn_pred_ibp <- vector(mode="list", length = length(Ts))
# names(list_kmn_pred_ibp) <- paste("T", Ts, sep = ".")
# list_kmn_pred_sp <- vector(mode="list", length = length(Ts))
# names(list_kmn_pred_sp) <- paste("T", Ts, sep = ".")
# 
# Kn = ncol(data_mat[,colSums(data_mat) > 0])
# 
# for (j in 1:length(Ts)){
#   hor <- Ts[j]
#   
#   # Poisson
#   kmn_chain_poiss <- generate_Kmn_chain_poiss(lambda_chain_poiss, alpha_chain_poiss,
#                                               theta_chain_poiss, M = hor, n = L)
#   
#   est_ci_pred_poiss <- matrix(NA, nrow = hor, ncol = 3)
#   # first column = lower bound
#   # second columns = medians
#   # third columns = upper bound
#   for (m in 1:hor){
#     est_ci_pred_poiss[m,] <- quantile(kmn_chain_poiss[m,], probs = c(0.025,0.5,0.975))
#   }
#   est_ci_pred_poiss <- list("medians" = est_ci_pred_poiss[,2],
#                             "lbs" = est_ci_pred_poiss[,1],
#                             "ubs" = est_ci_pred_poiss[,3])
#   
#   list_kmn_pred_poiss[[paste0("T.",hor)]] <- est_ci_pred_poiss
#   
#   # Negative Binomial
#   kmn_chain_negbin <- generate_Kmn_chain_negbin(nstar_chain_negbin, p_chain_negbin,
#                                                 alpha_chain_negbin, theta_chain_negbin,
#                                                 M = hor, n = L, Kn)
#   
#   est_ci_pred_negbin <- matrix(NA, nrow = hor, ncol = 3)
#   # first column = lower bound
#   # second columns = medians
#   # third columns = upper bound
#   for (m in 1:hor){
#     est_ci_pred_negbin[m,] <- quantile(kmn_chain_negbin[m,], probs = c(0.025,0.5,0.975))
#   }
#   est_ci_pred_negbin <- list("medians" = est_ci_pred_negbin[,2],
#                              "lbs" = est_ci_pred_negbin[,1],
#                              "ubs" = est_ci_pred_negbin[,3])
#   
#   list_kmn_pred_negbin[[paste0("T.",hor)]] <- est_ci_pred_negbin
#   
#   
#   # IBP + Gamma
#   kmn_chain_ibp <- generate_Kmn_chain_gamma_ibp(a_chain_ibp, b_chain_ibp, alpha_chain_ibp,
#                                                 theta_chain_ibp, M = hor, n = L, Kn)
#   
#   est_ci_pred_ibp <- matrix(NA, nrow = hor, ncol = 3)
#   # first column = lower bound
#   # second columns = medians
#   # third columns = upper bound
#   for (m in 1:hor){
#     est_ci_pred_ibp[m,] <- quantile(kmn_chain_ibp[m,], probs = c(0.025,0.5,0.975))
#   }
#   est_ci_pred_ibp <- list("medians" = est_ci_pred_ibp[,2],
#                           "lbs" = est_ci_pred_ibp[,1],
#                           "ubs" = est_ci_pred_ibp[,3])
#   
#   list_kmn_pred_ibp[[paste0("T.",hor)]] <- est_ci_pred_ibp
#   
#   # SB-SP
#   kmn_chain_sp <- generate_Kmn_chain_gamma_ibp(a_chain = c_chain_sp + 1,
#                                                b_chain = beta_chain_sp*(1-alpha_chain_sp)/alpha_chain_sp, 
#                                                alpha_chain = alpha_chain_sp,
#                                                theta_chain = 1 - alpha_chain_sp,
#                                                M = hor, n = L, Kn)
#   
#   est_ci_pred_sp <- matrix(NA, nrow = hor, ncol = 3)
#   # first column = lower bound
#   # second columns = medians
#   # third columns = upper bound
#   for (m in 1:hor){
#     est_ci_pred_sp[m,] <- quantile(kmn_chain_sp[m,], probs = c(0.025,0.5,0.975))
#   }
#   est_ci_pred_sp <- list("medians" = est_ci_pred_sp[,2],
#                          "lbs" = est_ci_pred_sp[,1],
#                          "ubs" = est_ci_pred_sp[,3])
#   
#   list_kmn_pred_sp[[paste0("T.",hor)]] <- est_ci_pred_sp
#   
#   
# }
# ######### Save results: CI for extrapolation (Poiss/NB/Gamma/SP)
# # Poisson
# saveRDS(list_kmn_pred_poiss, "bci_results/bci_op_small_ci_poiss.rds")
# # Negative Binomial
# saveRDS(list_kmn_pred_negbin, "bci_results/bci_op_small_ci_negbin.rds")
# # Gamma IBP
# saveRDS(list_kmn_pred_ibp, "bci_results/bci_op_small_ci_ibp.rds")
# # SB-SP
# saveRDS(list_kmn_pred_sp, "bci_results/bci_op_small_ci_sp.rds")


###### 7.d) Prediction on different horizons (results) - small-set new subjects ########
list_kmn_pred_poiss <- readRDS(file = "bci_results/bci_op_small_ci_poiss.rds")
list_kmn_pred_negbin <- readRDS(file = "bci_results/bci_op_small_ci_negbin.rds")
list_kmn_pred_ibp <- readRDS(file = "bci_results/bci_op_small_ci_ibp.rds")
list_kmn_pred_sp <- readRDS(file = "bci_results/bci_op_small_ci_sp.rds")

gg_kmn_pred_all  <- vector(mode="list", length = length(Ts))

for (j in 1:length(Ts)){
  hor <- Ts[j]
  
  gg_kmn_pred_all[[j]] <- plot_Kmn_median_pred_all(
    data_list = data_list,
    ci_poiss = list_kmn_pred_poiss[[paste0("T.",hor)]], 
    ci_negbin = list_kmn_pred_negbin[[paste0("T.",hor)]],
    ci_ibp = list_kmn_pred_ibp[[paste0("T.",hor)]],
    n_avg = 100) +
    ggtitle(paste0("# new observations = ", hor)) + theme(plot.title = element_text(size=12))
  
  
  if (j != 1){
    gg_kmn_pred_all[[j]] <- gg_kmn_pred_all[[j]] +
      theme(#axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        axis.title.y = element_blank() )
  }
  
}

# Print plots
fig_bci_pred_ <- wrap_plots(gg_kmn_pred_all, nrow = 1, ncol = length(Ts)) +   plot_layout(guides = "collect") & theme(legend.position = 'right') & xlab(NULL) & theme(plot.margin = margin(5.5, 5.5, 0, 5.5))

fig_bci_pred <- wrap_elements(panel = fig_bci_pred_) +
  labs(tag = "# observations") +
  theme(
    plot.tag = element_text(size = rel(1)),
    plot.tag.position = "bottom") 

ggsave(filename = "Plots_paper/plot_bci_op_small_pred.png", width = 8, height = 5, dpi = 300, units = "in", device='png')



###### 8) Richness
gg_ntilde_poiss_long <- gather(gg_ntilde_poiss, training, estimate, 
                               paste0("N.",Ns[1]):paste0("N.",Ns[length(Ns)]), 
                               factor_key=TRUE) %>%
  add_column(Model= "BBmixP", N = rep(Ns, each = nrow(gg_ntilde_poiss)))


gg_ntilde_negbin_long <- gather(gg_ntilde_negbin, training, estimate, 
                                paste0("N.",Ns[1]):paste0("N.",Ns[length(Ns)]), 
                                factor_key=TRUE) %>%
  add_column(Model= "BBmixNB", N = rep(Ns, each = nrow(gg_ntilde_negbin)))


joint_total_long <- bind_rows(gg_ntilde_poiss_long, gg_ntilde_negbin_long) %>%
  mutate(Model = fct_relevel(Model, c("BBmixP", "BBmixNB"))) 

# plot
ggplot(joint_total_long, aes(x = estimate, color = Model)) + 
  #geom_density(alpha = 0.5, linewidth = 0.8) +
  stat_density(aes(x=estimate, colour=Model),
               geom="line",position="identity") +
  geom_vline(xintercept = ncol(data_mat), color="grey", linetype="dashed", linewidth=0.8) + # observed features in the whole dataset
  facet_wrap(~N, labeller = labeller(N = ~ paste("n = ", .x)), scales = "free", nrow = 1) +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 12))+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(breaks = pretty_breaks()) +
  xlab("# distinct features") + rremove("ylab") +
  scale_color_manual(values = c("BBmixP" = "forestgreen",
                                "BBmixNB" = "royalblue1")) 

ggsave(filename = "Plots_paper/plot_bci_richness.png", width = 10, height = 4, dpi = 300, units = "in", device='png')
