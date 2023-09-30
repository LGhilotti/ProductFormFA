####
##### MODEL 1: the homogeneous model ############################
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

##### Single dataset -> Ntilde (Poiss/NB) and extrapolation (Poiss/NB/Gamma) #####

###### 1) Read results:  MCMC convergence ####################
load(file = "chiu_model_simulation/m1/m1_params_poiss.Rda")
load(file =  "chiu_model_simulation/m1/m1_params_negbin.Rda")
#load(file =  "chiu_model_simulation/m1/m1_params_negbin_prior.Rda")
load(file =  "chiu_model_simulation/m1/m1_params_ibp.Rda" )
#load(file =  "chiu_model_simulation/m1/m1_params_sp.Rda" )

###### 2) Read results: samples from limiting distributions (Poiss/NB) ##############
load(file = "chiu_model_simulation/m1/m1_ntilde_poiss.Rda")
load(file = "chiu_model_simulation/m1/m1_ntilde_negbin.Rda")
#load(file = "chiu_model_simulation/m1/m1_ntilde_negbin_prior.Rda")

list_kmn_pred_test_poiss <- readRDS(file = "chiu_model_simulation/m1/m1_ci_poiss.rds")

###### 3) Read the data ###############################
data_mat <- readRDS(file = "chiu_model_simulation/m1/m1_data_mat.rds")
data_list <- create_features_list(data_mat)
L <- nrow(data_mat)
Ms <- unique(sapply(list_kmn_pred_test_poiss, function(l) length(l$medians)))
Ns <- L - Ms

Nbars <- c(200,400,600)
c_fr <- 10

###### 4) Check MCMC convergence###################

# Poisson 

###### check mcmc mixing poisson Nbar = emp
params_poiss_ <- params_poiss[c(paste0("alpha:N.",Ns,":Nbar.emp"), 
                                paste0("theta:N.",Ns,":Nbar.emp")) ]
samples_poiss <- mcmc.list(mcmc(params_poiss_))
samples_ggs_poiss <- ggs(samples_poiss, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_poiss) + 
  facet_wrap(~Parameter, nrow = length(Ns), scales = "free")

effectiveSize(params_poiss_)

###### check mcmc mixing poisson Nbar = 200
params_poiss_ <- params_poiss[c(paste0("alpha:N.",Ns,":Nbar.200"), 
                                 paste0("theta:N.",Ns,":Nbar.200")) ]
samples_poiss <- mcmc.list(mcmc(params_poiss_))
samples_ggs_poiss <- ggs(samples_poiss, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_poiss) + 
  facet_wrap(~Parameter, nrow = length(Ns), scales = "free")

effectiveSize(params_poiss_)

###### check mcmc mixing poisson Nbar = 400
params_poiss_ <- params_poiss[c(paste0("alpha:N.",Ns,":Nbar.400"), 
                                paste0("theta:N.",Ns,":Nbar.400")) ]
samples_poiss <- mcmc.list(mcmc(params_poiss_))
samples_ggs_poiss <- ggs(samples_poiss, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_poiss) + 
  facet_wrap(~Parameter, nrow = length(Ns), scales = "free")

effectiveSize(params_poiss_)

###### check mcmc mixing poisson Nbar = 600
params_poiss_ <- params_poiss[c(paste0("alpha:N.",Ns,":Nbar.600"), 
                                paste0("theta:N.",Ns,":Nbar.600")) ]
samples_poiss <- mcmc.list(mcmc(params_poiss_))
samples_ggs_poiss <- ggs(samples_poiss, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_poiss) + 
  facet_wrap(~Parameter, nrow = length(Ns), scales = "free")

effectiveSize(params_poiss_)

# Negative Binomial

###### check mcmc mixing negbin (fixed) Nbar = 200
params_negbin_ <- params_negbin[c(paste0("alpha:N.",Ns,":Nbar.200"), 
                                paste0("theta:N.",Ns,":Nbar.200")) ]
samples_negbin <- mcmc.list(mcmc(params_negbin_))
samples_ggs_negbin <- ggs(samples_negbin, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_negbin) + 
  facet_wrap(~Parameter, nrow = length(Ns), scales = "free")

effectiveSize(params_negbin_)

###### check mcmc mixing negbin (fixed) Nbar = 400
params_negbin_ <- params_negbin[c(paste0("alpha:N.",Ns,":Nbar.400"), 
                                  paste0("theta:N.",Ns,":Nbar.400")) ]
samples_negbin <- mcmc.list(mcmc(params_negbin_))
samples_ggs_negbin <- ggs(samples_negbin, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_negbin) + 
  facet_wrap(~Parameter, nrow = length(Ns), scales = "free")

effectiveSize(params_negbin_)

###### check mcmc mixing negbin (fixed) Nbar = 600
params_negbin_ <- params_negbin[c(paste0("alpha:N.",Ns,":Nbar.600"), 
                                  paste0("theta:N.",Ns,":Nbar.600")) ]
samples_negbin <- mcmc.list(mcmc(params_negbin_))
samples_ggs_negbin <- ggs(samples_negbin, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_negbin) + 
  facet_wrap(~Parameter, nrow = length(Ns), scales = "free")

effectiveSize(params_negbin_)

# Gamma ibp

###### check mcmc mixing Gamma ibp
params_ibp_ <- params_ibp[c(paste0("alpha:N.",Ns), 
                                  paste0("theta:N.",Ns)) ]
samples_ibp <- mcmc.list(mcmc(params_ibp_))
samples_ggs_ibp <- ggs(samples_ibp, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_ibp) + 
  facet_wrap(~Parameter, nrow = length(Ns), scales = "free")

effectiveSize(params_ibp_)



####### 5) Plot limiting distributions (Poiss/NB) ##############
labels_comb_bb <- paste(rep(paste("N", Ns, sep = "."), each = length(Nbars)+1),
                     c(paste("Nbar", Nbars, sep = "."),"Nbar.emp"), sep=":")


gg_ntilde_poiss_long <- gather(gg_ntilde_poiss, training, estimate, 
                               labels_comb_bb, 
                               factor_key=TRUE) %>%
  add_column(Model= "BBmixP", 
             N = rep(Ns, each = nrow(gg_ntilde_poiss)*(length(Nbars)+1)),
             Nbar = rep(rep(c(Nbars,"emp"), each = nrow(gg_ntilde_poiss)), length(Ns) ) )


gg_ntilde_negbin_long <- gather(gg_ntilde_negbin, training, estimate, 
                                labels_comb_bb, 
                                factor_key=TRUE) %>%
  add_column(Model= "BBmixNB", 
             N = rep(Ns, each = nrow(gg_ntilde_negbin)*(length(Nbars)+1)),
             Nbar = rep(rep(c(Nbars,"emp"), each = nrow(gg_ntilde_negbin)), length(Ns) ))


joint_total_long <- bind_rows(gg_ntilde_poiss_long, gg_ntilde_negbin_long) %>%
  mutate(Model = fct_relevel(Model, c("BBmixP", "BBmixNB"))) 

# plot
ggplot(joint_total_long, aes(x = estimate, color = Model)) + 
  #geom_density(alpha = 0.5, linewidth = 0.8) +
  stat_density(aes(x=estimate, colour=Model),
               geom="line",position="identity") +
  geom_vline(xintercept = 500, color="black", linetype="dashed", linewidth=0.8) +
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

ggsave(filename = "Plots_paper/plot_m1_richness.png", width = 10, height = 4, dpi = 300, units = "in", device='png')


######## 6) Plot Extrapolation curve (Poiss/NB/Gamma) ###############
list_kmn_pred_test_poiss <- readRDS(file = "chiu_model_simulation/m1/m1_ci_poiss.rds")
list_kmn_pred_test_negbin <- readRDS(file = "chiu_model_simulation/m1/m1_ci_negbin.rds")
list_kmn_pred_test_ibp <- readRDS(file = "chiu_model_simulation/m1/m1_ci_ibp.rds")

# Poisson
list_kmn_pred_test_poiss <- lapply(list_kmn_pred_test_poiss, function(x) as_tibble(x))
list_kmn_pred_test_poiss <- lapply(1:length(list_kmn_pred_test_poiss), 
                                    function(x) (list_kmn_pred_test_poiss[[x]] %>%
                                                   add_column(Setting = names(list_kmn_pred_test_poiss[x]),
                                                              t = 1:nrow(list_kmn_pred_test_poiss[[x]])) %>%
                                                   extract(Setting, c("N","Nbar"), "N\\.([[:digit:]]+)\\:Nbar\\.([[:alnum:]]+)")
                                    ) )

df_pred_poiss <- bind_rows(list_kmn_pred_test_poiss) %>%
  add_column(Model = "BBmixP") 

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
  add_column(Model = "BBmixNB")

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
  add_column(Model = "3IBPmix",
             Nbar = "Not applicable")

df_pred_ibp$N <- as.integer(df_pred_ibp$N)
df_pred_ibp <- df_pred_ibp %>%
  mutate( t = t + N )

joint_df_pred_bayes <- rbind(df_pred_poiss,df_pred_negbin,df_pred_ibp)

## 6.b) Good-Toulmin prediction
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
# saveRDS(list_kmn_pred_test_gt, "chiu_model_simulation/m1/m1_gt_prediction.rds")

# read
list_kmn_pred_test_gt <- readRDS(file = "chiu_model_simulation/m1/m1_gt_prediction.rds")

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

list_chao_pred <- readRDS(file = "chiu_model_simulation/m1/m1_chao_rare.rds")
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
  filter(t <= L)


## 6.c) Observed sample <-> cts on the full sample
obs_sample <- sapply(2:L, function(n) ncol(data_mat[1:n,colSums(data_mat[1:n,]) > 0])   )
obs_sample <- data.frame(t = 0:L, 
                         obs = c(0, sum(data_mat[1,]) , obs_sample))

joint_df_pred_bayes <- joint_df_pred_bayes %>%
  mutate(medians = medians + obs_sample[N+1,]$obs,
         lbs = lbs + obs_sample[N+1,]$obs,
         ubs = ubs + obs_sample[N+1,]$obs)

###### 6.final) PLOT prediction ####
joint_df_pred_bayes_plot <- joint_df_pred_bayes %>%
  filter(Nbar %in% c("Not applicable", "emp"))
temp <- tibble(N = Ns, xvalues = Ns)

ggplot(joint_df_pred_bayes_plot, aes(t,medians, color = Model)) +
  geom_line(linetype = "dashed", linewidth = 0.8) +
  #geom_ribbon(aes(ymin = lbs, ymax = ubs), alpha = 0.1) +
  geom_line(data = obs_sample, aes(t, obs), color="black", linetype="solid", linewidth=0.4) +
  geom_line(data = df_pred_gt, aes(t, value), linetype = "dashed") +
  geom_line(data = df_pred_chao, aes(t, value), linetype = "dashed") +
  facet_wrap(.~ N, scales = "free_x") +
  geom_vline(data = temp, mapping =  aes(xintercept = xvalues) , linetype = "dashed", color = "grey") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_tableau() 

ggsave(filename = "Plots_paper/plot_m1_prediction.png", width = 10, height = 4, dpi = 300, units = "in", device='png')

##### 8) Richness: comparison with frequentist estimators ##########
table_richness <- joint_total_long %>% group_by(N, Nbar, Model) %>%
  summarise(estimator = mean(estimate), .groups = 'drop') %>%
  #summarise(estimator = quantile(estimate, probs = 0.25), .groups = 'drop') %>%
  mutate(Model_spec = paste0(Model,":Nbar.",Nbar) )

for (j in 1:length(Ns)){
  N <- Ns[j]
  M <- L - N
  
  train_mat <- data_mat[1:N,]
  
  # Beta-Binomial estimator
  table_richness <- table_richness %>%
    add_row(N = N, Nbar = NA, Model = "BetaBin", 
            estimator = beta_binomial_estimator(train_mat),
            Model_spec = "BetaBin")
  
  # Chao2 estimator
  table_richness <- table_richness %>%
    add_row(N = N, Nbar = NA, Model = "Chao2", 
            estimator = chao2_estimator(train_mat),
            Model_spec = "Chao2")
  
  # Improved Chao2 estimator
  table_richness <- table_richness %>%
    add_row(N = N, Nbar = NA, Model = "ImpChao2", 
            estimator = improved_chao2_estimator(train_mat),
            Model_spec = "ImpChao2")
  
  
}

table_richness <- arrange(table_richness,N,Nbar,Model)

table_err_richness <- table_richness %>% 
  mutate(err_estimate = abs(estimator - 500))

# plots estimator
table_richness$Nbar[is.na(table_richness$Nbar)] <- "Not applicable" 
ggplot(table_richness, aes( y=estimator, x=Model, shape = Nbar)) + 
  geom_point() + 
  facet_wrap(.~N, scales = "free_x", nrow = 1) +
  theme_light() + 
  geom_hline(aes(yintercept = 500), linetype = "dashed") +
  theme(legend.position = "top") +
  ylab("Estimate") +
  scale_y_continuous(breaks = pretty_breaks()) +
  rremove("xlab") 

ggsave(filename = "Plots_paper/plot_m1_richness.png", width = 10, height = 4, dpi = 300, units = "in", device='png')

# plots error
ggplot(table_err_richness, aes( y=err_estimate, x=Model, color=Model_spec)) + 
  geom_point(size = 3) + 
  facet_wrap(~N, labeller = labeller(N = ~ paste("n = ", .x)), scales = "free", nrow = 1) +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 12))+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  #theme(legend.position = "none") +
  scale_y_continuous(breaks = pretty_breaks()) +
  rremove("xlab") 



##### 9) Rarefaction curve and comparison with Chao bands #######

## 9.a) Bayesian estimates

list_kn_rarefaction_poiss <- readRDS(file = "chiu_model_simulation/m1/m1_ci_rare_poiss.rds")
list_kn_rarefaction_negbin <- readRDS(file = "chiu_model_simulation/m1/m1_ci_rare_negbin.rds")
list_kn_rarefaction_ibp <- readRDS(file = "chiu_model_simulation/m1/m1_ci_rare_ibp.rds")

# Poisson
list_kn_rarefaction_poiss <- lapply(list_kn_rarefaction_poiss, function(x) as_tibble(x))
list_kn_rarefaction_poiss <- lapply(1:length(list_kn_rarefaction_poiss), 
                                    function(x) (list_kn_rarefaction_poiss[[x]] %>%
                                                   add_column(Setting = names(list_kn_rarefaction_poiss[x]),
                                                              t = 1:nrow(list_kn_rarefaction_poiss[[x]])) %>%
                                                   extract(Setting, c("N","Nbar"), "N\\.([[:digit:]]+)\\:Nbar\\.([[:alnum:]]+)") 
                                                   ) )

df_rare_poiss <- bind_rows(list_kn_rarefaction_poiss) %>%
  add_column(Model = "BBmixP") 


# Negbin
list_kn_rarefaction_negbin <- lapply(list_kn_rarefaction_negbin, function(x) as_tibble(x))
list_kn_rarefaction_negbin <- lapply(1:length(list_kn_rarefaction_negbin), 
                                    function(x) (list_kn_rarefaction_negbin[[x]] %>%
                                                   add_column(Setting = names(list_kn_rarefaction_negbin[x]),
                                                              t = 1:nrow(list_kn_rarefaction_negbin[[x]])) %>%
                                                   extract(Setting, c("N","Nbar"), "N\\.([[:digit:]]+)\\:Nbar\\.([[:alnum:]]+)") 
                                    ) )

df_rare_negbin <- bind_rows(list_kn_rarefaction_negbin) %>%
  add_column(Model = "BBmixNB")

# Gamma IBP
list_kn_rarefaction_ibp <- lapply(list_kn_rarefaction_ibp, function(x) as_tibble(x))
list_kn_rarefaction_ibp <- lapply(1:length(list_kn_rarefaction_ibp), 
                                    function(x) (list_kn_rarefaction_ibp[[x]] %>%
                                                   add_column(Setting = names(list_kn_rarefaction_ibp[x]),
                                                              t = 1:nrow(list_kn_rarefaction_ibp[[x]])) %>%
                                                   extract(Setting, c("N"), "N\\.([[:digit:]]+)") 
                                    ) )

df_rare_ibp <- bind_rows(list_kn_rarefaction_ibp) %>%
  add_column(Model = "3IBPmix",
             Nbar = "Not applicable")

joint_df_rare_bayes <- rbind(df_rare_poiss,df_rare_negbin,df_rare_ibp)


## 9.b) Chao estimates

list_chao_rare <- readRDS(file = "chiu_model_simulation/m1/m1_chao_rare.rds")
list_chao_rare <- lapply(list_chao_rare, function(x) as_tibble(x))
list_chao_rare <- lapply(1:length(list_chao_rare), 
                         function(x) (list_chao_rare[[x]] %>%
                                        add_column(Setting = names(list_chao_rare[x])) %>%
                                        extract(Setting, c("N"), "N\\.([[:digit:]]+)") 
                                    ) )

df_rare_chao <- bind_rows(list_chao_rare) %>%
  add_column(Model = "Chao",
             Nbar = "Not applicable")

# Joint df with bayes and chao
joint_df_rare <- rbind(joint_df_rare_bayes, df_rare_chao)
joint_df_rare$N <- as.integer(joint_df_rare$N)

joint_df_rare_plot <- joint_df_rare  %>%
  filter((N +1)> t)

## 9.c) Observed sample
obs_rarefactions <- vector(mode="list", length = length(Ns))
names(obs_rarefactions) <- paste0("N.",Ns)
n_avg <- 100

for (j in 1:length(Ns)){
  N <- Ns[j]
  train_mat <- data_mat[1:N,]
  train_list <- create_features_list(train_mat)
  
  cum_nfeat <- rep(0, N)
  
  for (h in 1:n_avg){
    ord <- sample(1:N)
    cum_nfeat <- cum_nfeat + sapply(1:N, function(n) length(unique(unlist(train_list[ord][1:n]))))
  }
  
  cum_nfeat <- cum_nfeat/n_avg
  
  obs_rarefactions[[paste0("N.",N)]] <- data.frame(t = 0:N, obs = c(0,cum_nfeat), N = N)
}

obs_rarefactions <- bind_rows(obs_rarefactions)

#### 9.final) PLOT rarefaction ####
joint_df_rare_plot <- joint_df_rare_plot %>%
  filter(Nbar %in% c("Not applicable", "emp"))
temp <- tibble(N = Ns, xvalues = Ns)
  
ggplot(joint_df_rare_plot, aes(t,medians, color = Model)) +
  geom_line(linetype = "dashed", linewidth = 0.8) +
  #geom_ribbon(aes(ymin = lbs, ymax = ubs, fill = Model ), alpha = 0.1) +
  geom_line(data = obs_rarefactions, aes(t, obs), color="black", linetype="solid", linewidth=0.4) +
  facet_wrap(.~ N, scales = "free_x") +
  #geom_vline(data = temp, mapping =  aes(xintercept = xvalues) , linetype = "dashed", color = "grey") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_tableau()

ggsave(filename = "Plots_paper/plot_m1_rarefaction.png", width = 10, height = 4, dpi = 300, units = "in", device='png')






##### Accuracy on multiple datasets #####

###### 1) Read results: limit distribution estimates #####
load(file = "chiu_model_simulation/m1/m1_avg_ntilde_poiss.Rda")
load(file = "chiu_model_simulation/m1/m1_avg_ntilde_negbin.Rda")

###### 2) Read results: quantities on accuracy #####
load(file = "chiu_model_simulation/m1/m1_obs_train.Rda")
load(file = "chiu_model_simulation/m1/m1_obs_new.Rda")

load(file = "chiu_model_simulation/m1/m1_est_new_poiss.Rda")
load(file = "chiu_model_simulation/m1/m1_est_new_negbin.Rda")
load(file = "chiu_model_simulation/m1/m1_est_new_ibp.Rda")
load(file = "chiu_model_simulation/m1/m1_est_last_gt.Rda")
load(file = "chiu_model_simulation/m1/m1_est_last_chao.Rda")

D <- nrow(avg_ntilde_poiss)

###### 3) Plot limit distribution estimates #####
avg_ntilde_poiss_long <- avg_ntilde_poiss %>% 
  pivot_longer(everything(), names_to = "training", values_to = "estimate") %>%
  add_column(Model = "Poiss", N = rep(Ns[1:length(Ns)], each = D))

avg_ntilde_negbin_long <- avg_ntilde_negbin %>%
  pivot_longer(everything(), names_to = "training", values_to = "estimate") %>%
  add_column(Model = "Neg-Bin", N = rep(Ns[1:length(Ns)], each = D))

joint_ntilde_long <- bind_rows(avg_ntilde_poiss_long, avg_ntilde_negbin_long)

# plots
ggplot(joint_ntilde_long, aes(x=Model, y=estimate)) +
  geom_boxplot() +
  facet_wrap(~N) +
  geom_hline(aes(yintercept = 500), linetype = "dashed") +
  theme_light() + 
  #theme(legend.position = "top") +
  rremove("xlab") +
  ylab("Estimate") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_color_tableau()



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
  add_column(Model = "Poiss", N = rep(Ns, each = D))

acc_alt_negbin_long <- acc_alt_negbin %>%
  pivot_longer( everything(), names_to = "training", values_to = "Accuracy") %>%
  add_column(Model = "NB", N = rep(Ns, each = D))

acc_alt_ibp_long <- acc_alt_ibp %>%
  pivot_longer( everything(), names_to = "training", values_to = "Accuracy") %>%
  add_column(Model = "IBP", N = rep(Ns, each = D))

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
  facet_wrap(~N) +
  theme_light() +
  rremove("xlab") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_color_tableau()
