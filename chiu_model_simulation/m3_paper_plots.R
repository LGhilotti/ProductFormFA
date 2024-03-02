####
##### MODEL 3: the broken stick model ############################
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
load(file = "chiu_model_simulation/m3/m3_params_poiss.Rda")
load(file =  "chiu_model_simulation/m3/m3_params_negbin.Rda")
#load(file =  "chiu_model_simulation/m3/m3_params_negbin_prior.Rda")
load(file =  "chiu_model_simulation/m3/m3_params_ibp.Rda" )
#load(file =  "chiu_model_simulation/m3/m3_params_sp.Rda" )

###### 2) Read results: samples from limiting distributions (Poiss/NB) ##############
load(file = "chiu_model_simulation/m3/m3_ntilde_poiss.Rda")
load(file = "chiu_model_simulation/m3/m3_ntilde_negbin.Rda")
#load(file = "chiu_model_simulation/m3/m3_ntilde_negbin_prior.Rda")

list_kmn_pred_test_poiss <- readRDS(file = "chiu_model_simulation/m3/m3_ci_poiss.rds")

###### 3) Read the data ###############################
data_mat <- readRDS(file = "chiu_model_simulation/m3/m3_data_mat.rds")
data_list <- create_features_list(data_mat)
L <- nrow(data_mat)
Ms <- unique(sapply(list_kmn_pred_test_poiss, function(l) length(l$medians)))
Ns <- L - Ms

Nbars <- c(200,400,600)

####### 4.a) Richness: Point plot asymptotic mean (number of features) ##############
labels_comb_bb <- paste(rep(paste("N", Ns, sep = "."), each = length(Nbars)+1),
                        c(paste("Nbar", Nbars, sep = "."),"Nbar.emp"), sep=":")


gg_ntilde_poiss_long <- gather(gg_ntilde_poiss, training, estimate, 
                               labels_comb_bb, 
                               factor_key=TRUE) %>%
  add_column(Model= "Poisson", 
             N = rep(Ns, each = nrow(gg_ntilde_poiss)*(length(Nbars)+1)),
             Nbar = rep(rep(c(Nbars,"emp"), each = nrow(gg_ntilde_poiss)), length(Ns) ) )


gg_ntilde_negbin_long <- gather(gg_ntilde_negbin, training, estimate, 
                                labels_comb_bb, 
                                factor_key=TRUE) %>%
  add_column(Model= "NegBin", 
             N = rep(Ns, each = nrow(gg_ntilde_negbin)*(length(Nbars)+1)),
             Nbar = rep(rep(c(Nbars,"emp"), each = nrow(gg_ntilde_negbin)), length(Ns) ))


joint_total_long <- bind_rows(gg_ntilde_poiss_long, gg_ntilde_negbin_long) %>%
  mutate(Model = fct_relevel(Model, c( "NegBin", "Poisson"))) 


table_richness <- joint_total_long %>% group_by(N, Nbar, Model) %>%
  summarise(estimator = mean(estimate),
            lb = quantile(estimate, probs = 0.025),
            ub = quantile(estimate, probs = 0.975), .groups = 'drop') %>%
  mutate(Model_spec = paste0(Model,":Nbar.",Nbar) )


# plots estimator
ggplot(table_richness, aes( y=estimator, x=Model, shape = Nbar)) + 
  geom_point() + 
  facet_wrap(.~N, scales = "free_x", nrow = 1) +
  theme_light() + 
  geom_hline(aes(yintercept = 500), linetype = "dashed") +
  theme(legend.position = "top") +
  ylab("Posterior mean of N") +
  scale_y_continuous(breaks = pretty_breaks()) +
  rremove("xlab") +
  scale_shape_discrete(name = "Prior mean of N") +
  theme(aspect.ratio = 1)

ggsave(filename = "Plots_paper/plot_m3_richness_point.pdf", width = 10, height = 4, dpi = 300, units = "in", device='pdf')


###### 4.b) Richness: asymptotic distributions for Empirical ##############
joint_emp_long <- joint_total_long %>%
  filter(Nbar == "emp") %>%
  select(-training)

ggplot(joint_emp_long, aes(x = estimate, color = Model)) +
  #geom_density(alpha = 0.5, linewidth = 0.8) +
  stat_density(aes(x=estimate, colour=Model),
               geom="line",position="identity") +
  geom_vline(aes(xintercept = 500), linetype="dashed") +
  facet_wrap(.~N, scales = "free_x") +
  theme_light() +
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  xlab("# distinct features") + rremove("ylab") +
  scale_color_tableau() +
  theme(aspect.ratio = 1)


ggsave(filename = "Plots_paper/plot_m3_richness_distribution.pdf", width = 10, height = 4, dpi = 300, units = "in", device='pdf')


######## 5) Extrapolation - Poisson, NegBin, Chao and GT ###############
list_kmn_pred_test_poiss <- readRDS(file = "chiu_model_simulation/m3/m3_ci_poiss.rds")
list_kmn_pred_test_negbin <- readRDS(file = "chiu_model_simulation/m3/m3_ci_negbin.rds")
list_kmn_pred_test_ibp <- readRDS(file = "chiu_model_simulation/m3/m3_ci_ibp.rds")

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

joint_df_pred_bayes <- rbind(df_pred_poiss,df_pred_negbin,df_pred_ibp)

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
# saveRDS(list_kmn_pred_test_gt, "chiu_model_simulation/m3/m3_gt_prediction.rds")

# read
list_kmn_pred_test_gt <- readRDS(file = "chiu_model_simulation/m3/m3_gt_prediction.rds")

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


## 6.c) Chao estimates (for prediction)

list_chao_pred <- readRDS(file = "chiu_model_simulation/m3/m3_chao_rare.rds")
list_chao_pred <- lapply(list_chao_pred, function(x) as_tibble(x))
list_chao_pred <- lapply(1:length(list_chao_pred), 
                         function(x) (list_chao_pred[[x]] %>%
                                        add_column(Setting = names(list_chao_pred[x])) %>%
                                        extract(Setting, c("N"), "N\\.([[:digit:]]+)") 
                         ) )

df_pred_chao <- bind_rows(list_chao_pred) %>%
  select(-c(lbs,ubs)) %>%
  rename(value = medians) %>%
  add_column(Model = "Chao",
             Nbar = "Not applicable")


df_pred_chao$N <- as.integer(df_pred_chao$N)

df_pred_chao <- df_pred_chao %>%
  filter(t <= L, t > N)


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
                                         levels = c("Gamma","Poisson","NegBin", "Chao", "GT"))
temp <- tibble(N = Ns, xvalues = Ns)

ggplot(joint_df_pred_bayes_plot, aes(x = t, y = medians, color = Model)) +
  geom_line() +
  #geom_ribbon(aes(ymin = lbs, ymax = ubs), alpha = 0.1) +
  geom_line(data = obs_sample, aes(t, obs), color="black", linetype="solid") +
  geom_line(data = df_pred_gt, aes(t, value)) +
  geom_line(data = df_pred_chao, aes(t, value)) +
  facet_wrap(.~ N, scales = "free_x") +
  geom_vline(data = temp, mapping =  aes(xintercept = xvalues) , linetype = "dashed", color = "grey") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_tableau() +
  theme(aspect.ratio = 1)

ggsave(filename = "Plots_paper/plot_m3_prediction.pdf", width = 10, height = 4, dpi = 300, units = "in", device='pdf')








##### Accuracy on multiple datasets #####

###### 1) Read results: limit distribution estimates #####
load(file = "chiu_model_simulation/m3/m3_avg_ntilde_poiss.Rda")
load(file = "chiu_model_simulation/m3/m3_avg_ntilde_negbin.Rda")

###### 2) Read results: quantities on accuracy #####
load(file = "chiu_model_simulation/m3/m3_obs_train.Rda")
load(file = "chiu_model_simulation/m3/m3_obs_new.Rda")

load(file = "chiu_model_simulation/m3/m3_est_new_poiss.Rda")
load(file = "chiu_model_simulation/m3/m3_est_new_negbin.Rda")
load(file = "chiu_model_simulation/m3/m3_est_new_ibp.Rda")
load(file = "chiu_model_simulation/m3/m3_est_last_gt.Rda")
load(file = "chiu_model_simulation/m3/m3_est_last_chao.Rda")

D <- nrow(avg_ntilde_poiss)

###### 3) Plot limit distribution estimates #####
avg_ntilde_poiss_long <- avg_ntilde_poiss %>% 
  pivot_longer(everything(), names_to = "training", values_to = "estimate") %>%
  add_column(Model = "Poisson", N = rep(Ns[1:length(Ns)], each = D))

avg_ntilde_negbin_long <- avg_ntilde_negbin %>%
  pivot_longer(everything(), names_to = "training", values_to = "estimate") %>%
  add_column(Model = "NegBin", N = rep(Ns[1:length(Ns)], each = D))

joint_ntilde_long <- bind_rows(avg_ntilde_poiss_long, avg_ntilde_negbin_long)

# plots
ggplot(joint_ntilde_long, aes(x=Model, y=estimate)) +
  geom_boxplot() +
  facet_wrap(.~N, scales = "free_x", nrow = 1) +
  geom_hline(aes(yintercept = 500), linetype = "dashed") +
  theme_light() + 
  #theme(legend.position = "top") +
  rremove("xlab") +
  ylab("# distinct features") +
  scale_y_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1)

ggsave(filename = "Plots_paper/plot_m3_richness_rep.pdf", width = 10, height = 4, dpi = 300, units = "in", device='pdf')



###### 4) Plot boxplots on accuracy #####
# define the dataframe with estimated number of new features
colnames(est_new_poiss) <- paste("N",Ns,sep = ".")
colnames(est_new_negbin) <- paste("N",Ns,sep = ".")
colnames(est_new_ibp) <- paste("N",Ns,sep = ".")

est_new_gt <- est_last_gt - obs_train
est_new_chao <- est_last_chao - obs_train

# compute accuracy
acc_alt_poiss <- compute_accuracy(obs_new, est_new_poiss, obs_train)

acc_alt_negbin <- compute_accuracy(obs_new, est_new_negbin, obs_train)

acc_alt_ibp <- compute_accuracy(obs_new, est_new_ibp, obs_train)

acc_alt_gt <- compute_accuracy(obs_new, est_new_gt, obs_train)

acc_alt_chao <- compute_accuracy(obs_new, est_new_chao, obs_train)

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

acc_alt_gt_long <- acc_alt_gt %>%
  pivot_longer( everything(), names_to = "training", values_to = "Accuracy") %>%
  add_column(Model = "GT", N = rep(Ns, each = D))

acc_alt_chao_long <- acc_alt_chao %>%
  pivot_longer( everything(), names_to = "training", values_to = "Accuracy") %>%
  add_column(Model = "Chao", N = rep(Ns, each = D))


joint_alt_long <- bind_rows(acc_alt_poiss_long, acc_alt_negbin_long,
                            acc_alt_ibp_long, acc_alt_gt_long, acc_alt_chao_long)

# plots
ggplot(joint_alt_long, aes(x = Model, y=Accuracy)) +
  geom_boxplot() +
  facet_wrap(.~N, scales = "free_x", nrow = 1) +
  theme_light() +
  rremove("xlab") +
  ylab("Error index") +
  scale_y_continuous(breaks = pretty_breaks())+
  theme(aspect.ratio = 1)

ggsave(filename = "Plots_paper/plot_m3_accuracy_rep.pdf", width = 10, height = 4, dpi = 300, units = "in", device='pdf')

