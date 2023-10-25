####
##### Unbounded-features scenario_ Polynomial (exponent: 1) ############################
####

rm(list=ls())
library(ggmcmc)
library(coda)
library(ggpubr)
library(scales)
library(tidyverse)
library(gridExtra)
library(grid)
library(patchwork)
library(ProductFormFA)

##### Single dataset -> Ntilde (Poiss/NB) and extrapolation (Poiss/NB/Gamma) #####

###### 1) Read results:  MCMC convergence ####################
load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_params_poiss.Rda")
load(file =  "unbounded_features_simulation/unb_poly_1/unb_poly_1_params_negbin.Rda")
load(file =  "unbounded_features_simulation/unb_poly_1/unb_poly_1_params_ibp.Rda" )
load(file =  "unbounded_features_simulation/unb_poly_1/unb_poly_1_params_sp.Rda" )

###### 2) Read results: samples from limiting distributions (Poiss/NB) ##############
load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_ntilde_poiss.Rda")
load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_ntilde_negbin.Rda")

###### 3) Read results: CI for extrapolation (Poiss/NB/Gamma) ################
list_kmn_pred_test_poiss <- readRDS(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_ci_poiss.rds")
list_kmn_pred_test_negbin <- readRDS(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_ci_negbin.rds")
list_kmn_pred_test_ibp <- readRDS(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_ci_ibp.rds")
list_kmn_pred_test_sp <- readRDS(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_ci_sp.rds")



###### 4) Read the data ###############################
data_mat <- readRDS(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_data_mat.rds")
data_list <- create_features_list(data_mat)
L <- nrow(data_mat)
Ms <- unique(sapply(list_kmn_pred_test_poiss, function(l) length(l$medians)))
Ns <- L - Ms


###### 5) Check MCMC convergence###################

# Poisson 

###### check mcmc mixing poisson Nbar=emp
params_poiss_ <- params_poiss[c(paste0("alpha:N.",Ns,":Nbar.emp" ), 
                                paste0("theta:N.",Ns,":Nbar.emp")) ]
samples_poiss <- mcmc.list(mcmc(params_poiss_))
samples_ggs_poiss <- ggs(samples_poiss, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_poiss) + 
  facet_wrap(~Parameter, nrow = length(Ns), scales = "free")

effectiveSize(params_poiss_)



# Negative Binomial

###### check mcmc mixing negbin (fixed) Nbar = emp
params_negbin_ <- params_negbin[c(paste0("alpha:N.",Ns,":Nbar.emp"), 
                                  paste0("theta:N.",Ns,":Nbar.emp")) ]
samples_negbin <- mcmc.list(mcmc(params_negbin_))
samples_ggs_negbin <- ggs(samples_negbin, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_negbin) + 
  facet_wrap(~Parameter, nrow = length(Ns), scales = "free")

effectiveSize(params_negbin_)



# Gamma ibp

###### check mcmc mixing Gamma ibp Nbar=emp
params_ibp_ <- params_ibp[c(paste0("alpha:N.",Ns), 
                            paste0("theta:N.",Ns)) ]
samples_ibp <- mcmc.list(mcmc(params_ibp_))
samples_ggs_ibp <- ggs(samples_ibp, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_ibp) + 
  facet_wrap(~Parameter, nrow = length(Ns), scales = "free")

effectiveSize(params_ibp_)


# SB-SP

###### check mcmc mixing SB_SP Nbar=emp
params_sp_ <- params_sp[c(paste0("alpha:N.",Ns), 
                          paste0("c:N.",Ns),
                          paste0("beta:N.",Ns)) ]
samples_sp <- mcmc.list(mcmc(params_sp_))
samples_ggs_sp <- ggs(samples_sp, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_sp) + 
  facet_wrap(~Parameter, nrow = length(Ns), scales = "free")

effectiveSize(params_sp_)



####### 6) Plot limiting distributions (Poiss/NB) ##############
labels_comb <- paste(rep(paste("N", Ns, sep = "."), each = length(Nbars)),
                     paste("Nbar", Nbars, sep = "."), sep=":")


gg_ntilde_poiss_long <- gather(gg_ntilde_poiss, training, estimate, 
                               labels_comb, 
                               factor_key=TRUE) %>%
  add_column(Model= "BBmixP", 
             N = rep(Ns, each = nrow(gg_ntilde_poiss)*length(Nbars)),
             Nbar = rep(rep(Nbars, each = nrow(gg_ntilde_poiss)), length(Ns) ) )


gg_ntilde_negbin_long <- gather(gg_ntilde_negbin, training, estimate, 
                                labels_comb, 
                                factor_key=TRUE) %>%
  add_column(Model= "BBmixNB", 
             N = rep(Ns, each = nrow(gg_ntilde_negbin)*length(Nbars)),
             Nbar = rep(rep(Nbars, each = nrow(gg_ntilde_negbin)), length(Ns) ))


joint_total_long <- bind_rows(gg_ntilde_poiss_long, gg_ntilde_negbin_long) %>%
  mutate(Model = fct_relevel(Model, c("BBmixP", "BBmixNB"))) 

# plot
ggplot(joint_total_long, aes(x = estimate, color = Model)) + 
  #geom_density(alpha = 0.5, linewidth = 0.8) +
  stat_density(aes(x=estimate, colour=Model),
               geom="line",position="identity") +
  #geom_vline(xintercept = 500, color="black", linetype="dashed", linewidth=0.8) +
  facet_wrap(~N+ Nbar, labeller = labeller(N = ~ paste("n = ", .x), Nbar = ~ paste("Nbar = ", .x)), scales = "free") +
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

ggsave(filename = "Plots_paper/plot_unb_poly_1_richness.png", width = 10, height = 4, dpi = 300, units = "in", device='png')


######## 7) Plot Extrapolation curve (Poiss/NB/Gamma) ################

gg_kmn_pred_test_all  <- vector(mode="list", length = length(Ns)*length(Nbars))
names(gg_kmn_pred_test_all) <- labels_comb

for (j in 1:length(Ns)){
  N <- Ns[j]
  M <- L - N
  
  train_mat <- data_mat[1:N,]
  test_mat <- data_mat[(N+1):L, ]
  # convert the binary matrix into list of features
  train_list <- create_features_list(train_mat)
  test_list <- create_features_list(test_mat)
  
  for (v in 1:length(Nbars)){
    Nbar <- Nbars[v]
    
    lab_comb <- paste0("N.",N,":Nbar.",Nbar)
    
    gg_kmn_pred_test_all[[lab_comb]] <- plot_Kmn_median_pred_and_test_all(
      train_list = train_list,
      test_list = test_list,
      ci_poiss = list_kmn_pred_test_poiss[[lab_comb]], 
      ci_negbin = list_kmn_pred_test_negbin[[lab_comb]],
      #ci_negbin_prior = list_kmn_pred_test_negbin_prior[[lab_comb]],
      ci_ibp = list_kmn_pred_test_ibp[[lab_comb]],
      n_avg = 100) +
      ggtitle(lab_comb) + theme(plot.title = element_text(size=12))
    
    
  }
  
  
  # if (j != 1){
  #   gg_kmn_pred_test_all[[j]] <- gg_kmn_pred_test_all[[j]] +
  #     theme(#axis.text.y = element_blank(),
  #       #axis.ticks.y = element_blank(),
  #       axis.title.y = element_blank() )
  # }
  
}

# Print plots
fig_unb_poly_1_pred_ <- wrap_plots(gg_kmn_pred_test_all, nrow = 3, ncol = 3) +   plot_layout(guides = "collect") & theme(legend.position = 'right') & xlab(NULL) & theme(plot.margin = margin(5.5, 5.5, 0, 5.5))

fig_unb_poly_1_pred <- wrap_elements(panel = fig_unb_poly_1_pred_) +
  labs(tag = "# observations") +
  theme(
    plot.tag = element_text(size = rel(1)),
    plot.tag.position = "bottom") 

ggsave(filename = "Plots_paper/plot_unb_poly_1_pred.png", width = 10, height = 4, dpi = 300, units = "in", device='png')





##### Accuracy on multiple datasets #####

###### 1) Read results: limit distribution estimates #####
load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_avg_ntilde_poiss.Rda")
load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_avg_ntilde_negbin.Rda")

###### 2) Read results: quantities on accuracy #####
load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_obs_train.Rda")
load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_obs_new.Rda")

load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_est_new_poiss.Rda")
load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_est_new_negbin.Rda")
load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_est_new_ibp.Rda")
load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_est_new_sp.Rda")

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



###### 4.1) Plot boxplots on accuracy (scaled) #####
ratio_scaled_poiss <- as.matrix(abs(obs_new - est_new_poiss)/obs_new)
acc_scaled_poiss <- as.data.frame(1 - pmin(ratio_scaled_poiss,1))

ratio_scaled_negbin <- as.matrix(abs(obs_new - est_new_negbin)/obs_new)
acc_scaled_negbin <- as.data.frame(1 - pmin(ratio_scaled_negbin,1))

ratio_scaled_ibp <- as.matrix(abs(obs_new - est_new_ibp)/obs_new)
acc_scaled_ibp <- as.data.frame(1 - pmin(ratio_scaled_ibp,1))

acc_scaled_poiss_long <- gather(acc_scaled_poiss, training, accuracy, 
                                paste0("N.",Ns[1]):paste0("N.",Ns[length(Ns)]), factor_key=TRUE) %>%
  add_column(model = "Poiss", N = rep(Ns, each = D))

acc_scaled_negbin_long <- gather(acc_scaled_negbin, training , accuracy, 
                                 paste0("N.",Ns[1]):paste0("N.",Ns[length(Ns)]), factor_key=TRUE) %>%
  add_column(model = "Neg-Bin", N = rep(Ns, each = D))

acc_scaled_ibp_long <- gather(acc_scaled_ibp, training , accuracy, 
                              paste0("N.",Ns[1]):paste0("N.",Ns[length(Ns)]), factor_key=TRUE) %>%
  add_column(model = "Gamma IBP", N = rep(Ns, each = D))

joint_scaled_long <- bind_rows(acc_scaled_poiss_long, acc_scaled_negbin_long, acc_scaled_ibp_long)

# plots
ggplot(joint_scaled_long, aes( y=accuracy, fill=model)) + 
  geom_boxplot() + 
  facet_wrap(~N) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(paste0("Model 1: accuracy (scaled) on D = ", D)) +
  scale_y_continuous(breaks = pretty_breaks()) 


###### 4.2) Plot boxplots on accuracy (alt) #####
acc_alt_poiss <- 1/(1 + abs(obs_new - est_new_poiss))

acc_alt_negbin <- 1/(1 + abs(obs_new - est_new_negbin))

acc_alt_ibp <- 1/(1 + abs(obs_new - est_new_ibp))

acc_alt_sp <- 1/(1 + abs(obs_new - est_new_sp))


acc_alt_poiss_long <- gather(acc_alt_poiss, training, accuracy, 
                             paste0("N.",Ns[1]):paste0("N.",Ns[length(Ns)]), factor_key=TRUE) %>%
  add_column(model = "Poiss", N = rep(Ns, each = D))

acc_alt_negbin_long <- gather(acc_alt_negbin, training , accuracy, 
                              paste0("N.",Ns[1]):paste0("N.",Ns[length(Ns)]), factor_key=TRUE) %>%
  add_column(model = "Neg-Bin", N = rep(Ns, each = D))

acc_alt_ibp_long <- gather(acc_alt_ibp, training , accuracy, 
                           paste0("N.",Ns[1]):paste0("N.",Ns[length(Ns)]), factor_key=TRUE) %>%
  add_column(model = "Gamma IBP", N = rep(Ns, each = D))

acc_alt_sp_long <- gather(acc_alt_sp, training , accuracy, 
                          paste0("N.",Ns[1]):paste0("N.",Ns[length(Ns)]), factor_key=TRUE) %>%
  add_column(model = "SB-SP", N = rep(Ns, each = D))

joint_alt_long <- bind_rows(acc_alt_poiss_long, acc_alt_negbin_long, 
                            acc_alt_ibp_long, acc_alt_sp_long)

# plots
ggplot(joint_alt_long, aes( y=accuracy, fill=model)) + 
  geom_boxplot() + 
  facet_wrap(~N) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle(paste0("Model 3: accuracy (alt) on D = ", D)) +
  scale_y_continuous(breaks = pretty_breaks()) 

