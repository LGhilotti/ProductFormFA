####
##### Unbounded with xi = 1 ############################
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
library(ggthemes)

##### Single dataset #####

###### 1) Read results:  MCMC convergence ####################
load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_params_poiss.Rda")
load(file =  "unbounded_features_simulation/unb_poly_1/unb_poly_1_params_negbin.Rda")
load(file =  "unbounded_features_simulation/unb_poly_1/unb_poly_1_params_ibp.Rda" )
load(file =  "unbounded_features_simulation/unb_poly_1/unb_poly_1_params_sp.Rda" )

###### 2) Read results: samples from limiting distributions (Poiss/NB) ##############
load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_ntilde_poiss.Rda")
load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_ntilde_negbin.Rda")

list_kmn_pred_test_poiss <- readRDS(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_ci_poiss.rds")

###### 3) Read the data ###############################
data_mat <- readRDS(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_data_mat.rds")
data_list <- create_features_list(data_mat)
L <- nrow(data_mat)
Ms <- unique(sapply(list_kmn_pred_test_poiss, function(l) length(l$medians)))
Ns <- L - Ms


###### 4) CHECK: MCMC convergence###################

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



######## 5) Extrapolation - Poisson, NegBin, Chao and GT ###############
list_kmn_pred_test_poiss <- readRDS(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_ci_poiss.rds")
list_kmn_pred_test_negbin <- readRDS(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_ci_negbin.rds")
list_kmn_pred_test_ibp <- readRDS(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_ci_ibp.rds")
list_kmn_pred_test_sp <- readRDS(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_ci_sp.rds")

# Poisson
list_kmn_pred_test_poiss <- lapply(list_kmn_pred_test_poiss, function(x) as_tibble(x))
list_kmn_pred_test_poiss <- lapply(1:length(list_kmn_pred_test_poiss), 
                                   function(x) (list_kmn_pred_test_poiss[[x]] %>%
                                                  add_column(Setting = names(list_kmn_pred_test_poiss[x]),
                                                             t = 1:nrow(list_kmn_pred_test_poiss[[x]])) %>%
                                                  extract(Setting, c("N","Nbar"), "N\\.([[:digit:]]+)\\:Nbar\\.([[:alnum:]]+)")
                                   ) )

df_pred_poiss <- bind_rows(list_kmn_pred_test_poiss) %>%
  add_column(Model = "Poisson") 

df_pred_poiss$N <- as.integer(df_pred_poiss$N)
df_pred_poiss <- df_pred_poiss %>%
  mutate( t = t + N )


# Negbin
list_kmn_pred_test_negbin <- lapply(list_kmn_pred_test_negbin, function(x) as_tibble(x))
list_kmn_pred_test_negbin <- lapply(1:length(list_kmn_pred_test_negbin), 
                                    function(x) (list_kmn_pred_test_negbin[[x]] %>%
                                                   add_column(Setting = names(list_kmn_pred_test_negbin[x]),
                                                              t = 1:nrow(list_kmn_pred_test_negbin[[x]])) %>%
                                                   extract(Setting, c("N","Nbar"), "N\\.([[:digit:]]+)\\:Nbar\\.([[:alnum:]]+)")
                                    ) )

df_pred_negbin <- bind_rows(list_kmn_pred_test_negbin) %>%
  add_column(Model = "NegBin")

df_pred_negbin$N <- as.integer(df_pred_negbin$N)
df_pred_negbin <- df_pred_negbin %>%
  mutate( t = t + N )

# Gamma IBP
list_kmn_pred_test_ibp <- lapply(list_kmn_pred_test_ibp, function(x) as_tibble(x))
list_kmn_pred_test_ibp <- lapply(1:length(list_kmn_pred_test_ibp), 
                                 function(x) (list_kmn_pred_test_ibp[[x]] %>%
                                                add_column(Setting = names(list_kmn_pred_test_ibp[x]),
                                                           t = 1:nrow(list_kmn_pred_test_ibp[[x]]) ) %>%
                                                extract(Setting, c("N"), "N\\.([[:digit:]]+)")
                                 ) )

df_pred_ibp <- bind_rows(list_kmn_pred_test_ibp) %>%
  add_column(Model = "Gamma",
             Nbar = "Not applicable")

df_pred_ibp$N <- as.integer(df_pred_ibp$N)
df_pred_ibp <- df_pred_ibp %>%
  mutate( t = t + N )

# SB-SP
list_kmn_pred_test_sp <- lapply(list_kmn_pred_test_sp, function(x) as_tibble(x))
list_kmn_pred_test_sp <- lapply(1:length(list_kmn_pred_test_sp), 
                                 function(x) (list_kmn_pred_test_sp[[x]] %>%
                                                add_column(Setting = names(list_kmn_pred_test_sp[x]),
                                                           t = 1:nrow(list_kmn_pred_test_sp[[x]]) ) %>%
                                                extract(Setting, c("N"), "N\\.([[:digit:]]+)")
                                 ) )

df_pred_sp <- bind_rows(list_kmn_pred_test_sp) %>%
  add_column(Model = "SB-Scaled",
             Nbar = "Not applicable")

df_pred_sp$N <- as.integer(df_pred_sp$N)
df_pred_sp <- df_pred_sp %>%
  mutate( t = t + N )

joint_df_pred_bayes <- rbind(df_pred_poiss,df_pred_negbin,df_pred_ibp, df_pred_sp)

## 5.b) Good-Toulmin prediction
# compute
# labels_comb_ibp <- paste("N", Ns, sep = ".")
# 
# list_kmn_pred_test_gt <- vector(mode="list", length = length(Ns))
# names(list_kmn_pred_test_gt) <- labels_comb_ibp
# 
# for (j in 1:length(Ns)){
#   
#   N <- Ns[j]
#   M <- L - N
#   
#   train_mat <- data_mat[1:N,]
#   # convert the binary matrix into list of features
#   train_list <- create_features_list(train_mat)
#   
#   lab_comb <- paste0("N.",N)
#   
#   # Compute SFS vector and CTS vector
#   sfs <- tabulate(colSums(train_mat))
#   
#   cts <- sapply(2:N, function(n) ncol(train_mat[1:n,colSums(train_mat[1:n,]) > 0])   )
#   cts <- c(0, sum(train_mat[1,]) , cts)
#   
#   list_kmn_pred_test_gt[[lab_comb]] <- predict_good_toulmin(N, M, sfs, cts, alternative = 0)$preds
#   
# }
# # Good-Toulmin predictions
# saveRDS(list_kmn_pred_test_gt, "unbounded_features_simulation/unb_poly_1/unb_poly_1_gt_prediction.rds")

# read
list_kmn_pred_test_gt <- readRDS(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_gt_prediction.rds")

list_kmn_pred_test_gt <- lapply(list_kmn_pred_test_gt, function(x) as_tibble(x))
list_kmn_pred_test_gt <- lapply(1:length(list_kmn_pred_test_gt), 
                                function(x) (list_kmn_pred_test_gt[[x]] %>%
                                               add_column(Setting = names(list_kmn_pred_test_gt[x]),
                                                          t = 0:(nrow(list_kmn_pred_test_gt[[x]]) -1)) %>%
                                               extract(Setting, c("N"), "N\\.([[:digit:]]+)") 
                                ) )

df_pred_gt <- bind_rows(list_kmn_pred_test_gt) %>%
  add_column(Model = "GT",
             Nbar = "Not applicable") 

df_pred_gt$N <- as.integer(df_pred_gt$N)

df_pred_gt <- df_pred_gt %>%
  filter(t > N)



## 5.c) Observed sample <-> cts on the full sample
obs_sample <- sapply(2:L, function(n) length(unique(unlist(data_list[1:n]))) )

obs_sample <- data.frame(t = 0:L, 
                         obs = c(0, sum(data_mat[1,]) , obs_sample))

joint_df_pred_bayes <- joint_df_pred_bayes %>%
  mutate(medians = medians + obs_sample[N+1,]$obs,
         lbs = lbs + obs_sample[N+1,]$obs,
         ubs = ubs + obs_sample[N+1,]$obs)


# plot
joint_df_pred_bayes_plot <- joint_df_pred_bayes %>%
  filter(Nbar %in% c("Not applicable", "emp"))
joint_df_pred_bayes_plot$Model <- factor(joint_df_pred_bayes_plot$Model, 
                                         levels = c("Gamma","Poisson","NegBin", "SB-Scaled", "GT"))
temp <- tibble(N = Ns, xvalues = Ns)

ggplot(joint_df_pred_bayes_plot, aes(x = t, y = medians, color = Model)) +
  geom_line(linetype = "dashed") +
  geom_ribbon(aes(ymin = lbs, ymax = ubs), alpha= 0.05) +
  geom_line(data = obs_sample, aes(t, obs), color="black", linetype="solid") +
  geom_line(data = df_pred_gt, aes(t, value)) +
  facet_wrap(.~ N, scales = "free_x") +
  geom_vline(data = temp, mapping =  aes(xintercept = xvalues) , linetype = "dashed", color = "grey") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_tableau() +
  theme(aspect.ratio = 1)

ggsave(filename = "Plots_paper/plot_unb_poly_1_prediction.pdf", width = 10, height = 4, dpi = 300, units = "in", device='pdf')








##### Error index on multiple datasets #####

###### 1) Read results: limit distribution estimates #####
load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_avg_ntilde_poiss.Rda")
load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_avg_ntilde_negbin.Rda")

###### 2) Read results: quantities on Error index #####
load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_obs_train.Rda")
load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_obs_new.Rda")

load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_est_new_poiss.Rda")
load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_est_new_negbin.Rda")
load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_est_new_ibp.Rda")
load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_est_new_sp.Rda")

load(file = "unbounded_features_simulation/unb_poly_1/unb_poly_1_est_last_gt.Rda")

D <- nrow(avg_ntilde_poiss)



###### 4) Plot boxplots on Error index #####
# define the dataframe with estimated number of new features
colnames(est_new_poiss) <- paste("N",Ns,sep = ".")
colnames(est_new_negbin) <- paste("N",Ns,sep = ".")
colnames(est_new_ibp) <- paste("N",Ns,sep = ".")
colnames(est_new_sp) <- paste("N",Ns,sep = ".")

est_new_gt <- est_last_gt - obs_train

# compute Error index
acc_alt_poiss <- compute_accuracy(obs_new, est_new_poiss, obs_train)

acc_alt_negbin <- compute_accuracy(obs_new, est_new_negbin, obs_train)

acc_alt_ibp <- compute_accuracy(obs_new, est_new_ibp, obs_train)

acc_alt_sp <- compute_accuracy(obs_new, est_new_sp, obs_train)

acc_alt_gt <- compute_accuracy(obs_new, est_new_gt, obs_train)


# handle the dataframes
acc_alt_poiss_long <- acc_alt_poiss %>%
  pivot_longer( everything(), names_to = "training", values_to = "Accuracy") %>%
  add_column(Model = "Poisson", N = rep(Ns, each = D))

acc_alt_negbin_long <- acc_alt_negbin %>%
  pivot_longer( everything(), names_to = "training", values_to = "Accuracy") %>%
  add_column(Model = "NegBin", N = rep(Ns, each = D))

acc_alt_ibp_long <- acc_alt_ibp %>%
  pivot_longer( everything(), names_to = "training", values_to = "Accuracy") %>%
  add_column(Model = "Gamma", N = rep(Ns, each = D))

acc_alt_sp_long <- acc_alt_sp %>%
  pivot_longer( everything(), names_to = "training", values_to = "Accuracy") %>%
  add_column(Model = "SB-Scaled", N = rep(Ns, each = D))

acc_alt_gt_long <- acc_alt_gt %>%
  pivot_longer( everything(), names_to = "training", values_to = "Accuracy") %>%
  add_column(Model = "GT", N = rep(Ns, each = D))



joint_alt_long <- bind_rows(acc_alt_poiss_long, acc_alt_negbin_long,
                            acc_alt_ibp_long,acc_alt_sp_long, acc_alt_gt_long)

# plots
ggplot(joint_alt_long, aes(x = Model, y=Accuracy)) +
  geom_boxplot() +
  facet_wrap(~N) +
  theme_light() +
  rremove("xlab") +
  ylab("Error index") +
  scale_y_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1)

ggsave(filename = "Plots_paper/plot_unb_poly_1_accuracy_rep.pdf", width = 10, height = 4, dpi = 300, units = "in", device='pdf')


###### 5) Comparison between Gamma and SB in terms of Error index #####
joint_alt_long_rest <- joint_alt_long %>%
  filter(Model %in% c("Gamma", "SB-Scaled"))

err_means <- joint_alt_long_rest %>% 
  group_by(Model) %>% 
  summarize(means=mean(Accuracy)) 

ggplot(joint_alt_long_rest, aes(x = Accuracy, color = Model)) +
  #geom_density(alpha = 0.5, linewidth = 0.8) +
  stat_density(aes(x=Accuracy, colour=Model),
               geom="line",position="identity") +
  geom_vline(data = err_means, aes(xintercept = means, color = Model), linetype="dashed") +
  facet_wrap(.~N, scales = "free_x") +
  theme_light() +
  scale_y_continuous(breaks = pretty_breaks()) +
  xlab("Error index") + rremove("ylab") +
  scale_color_tableau() 
