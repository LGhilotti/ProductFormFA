####
##### Plants ############################
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
load(file = "mazziotta2016_application/plants/mazz_plants_op_params_poiss.Rda")
load(file =  "mazziotta2016_application/plants/mazz_plants_op_params_negbin.Rda")
load(file =  "mazziotta2016_application/plants/mazz_plants_op_params_ibp.Rda" )

###### 2) Read results: samples from limiting distributions (Poiss/NB) ##############
load(file = "mazziotta2016_application/plants/mazz_plants_op_ntilde_poiss.Rda")
load(file = "mazziotta2016_application/plants/mazz_plants_op_ntilde_negbin.Rda")

list_kmn_pred_poiss <- readRDS(file = "mazziotta2016_application/plants/mazz_plants_op_ci_pred_poiss.rds")

###### 3) Read the data ###############################
data_mat <- readRDS(file = "mazziotta2016_application/plants/mazz_plants_data_mat.rds")
data_list <- create_features_list(data_mat)
L <- nrow(data_mat)


###### 4) CHECK: MCMC convergence###################

# Poisson 

###### check mcmc mixing poisson Nbar=emp
samples_poiss <- mcmc.list(mcmc(params_poiss))
samples_ggs_poiss <- ggs(samples_poiss, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_poiss) + 
  facet_wrap(~Parameter, scales = "free")

effectiveSize(params_poiss)



# Negative Binomial

###### check mcmc mixing negbin (fixed) Nbar = emp
samples_negbin <- mcmc.list(mcmc(params_negbin))
samples_ggs_negbin <- ggs(samples_negbin, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_negbin) + 
  facet_wrap(~Parameter,  scales = "free")

effectiveSize(params_negbin)



# Gamma ibp

###### check mcmc mixing Gamma ibp Nbar=emp
samples_ibp <- mcmc.list(mcmc(params_ibp))
samples_ggs_ibp <- ggs(samples_ibp, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_ibp) + 
  facet_wrap(~Parameter, scales = "free")

effectiveSize(params_ibp)


###### 4.b) Richness: asymptotic distributions for Empirical ##############

gg_ntilde_poiss_long <- gg_ntilde_poiss %>%
  add_column(Model= "Poisson" ) %>%
  rename(estimate = Nbar.emp)


gg_ntilde_negbin_long <- gg_ntilde_negbin %>%
  add_column(Model= "NegBin") %>%
  rename(estimate = Nbar.emp)


joint_emp_long <- bind_rows(gg_ntilde_poiss_long, gg_ntilde_negbin_long) %>%
  mutate(Model = fct_relevel(Model, c( "NegBin", "Poisson"))) 


ggplot(joint_emp_long, aes(x = estimate, color = Model)) +
  #geom_density(alpha = 0.5, linewidth = 0.8) +
  stat_density(aes(x=estimate, colour=Model),
               geom="line",position="identity") +
  theme_light() +
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  xlab("# distinct features") + rremove("ylab") +
  scale_color_tableau() +
  theme(aspect.ratio = 1)


ggsave(filename = "Plots_paper/plot_plants_op_richness_distribution.png", width = 10, height = 4, dpi = 300, units = "in", device='png')


######## 5) Extrapolation - Poisson, NegBin, Chao and GT ###############
list_kmn_pred_poiss <- readRDS(file = "mazziotta2016_application/plants/mazz_plants_op_ci_pred_poiss.rds")
list_kmn_pred_negbin <- readRDS(file = "mazziotta2016_application/plants/mazz_plants_op_ci_pred_negbin.rds")
list_kmn_pred_ibp <- readRDS(file = "mazziotta2016_application/plants/mazz_plants_op_ci_pred_ibp.rds")

# Poisson
list_kmn_pred_poiss <- lapply(list_kmn_pred_poiss, function(x) as_tibble(x))

df_pred_poiss <- bind_rows(list_kmn_pred_poiss) %>%
  add_column(N = L, t = 1:nrow(list_kmn_pred_poiss[[1]]),
             Model = "Poisson") 

df_pred_poiss$N <- as.integer(df_pred_poiss$N)
df_pred_poiss <- df_pred_poiss %>%
  mutate( t = t + N )


# Negbin
list_kmn_pred_negbin <- lapply(list_kmn_pred_negbin, function(x) as_tibble(x))

df_pred_negbin <- bind_rows(list_kmn_pred_negbin) %>%
  add_column(N = L, t = 1:nrow(list_kmn_pred_negbin[[1]]),
             Model = "NegBin")

df_pred_negbin$N <- as.integer(df_pred_negbin$N)
df_pred_negbin <- df_pred_negbin %>%
  mutate( t = t + N )

# Gamma IBP
list_kmn_pred_ibp <- lapply(list_kmn_pred_ibp, function(x) as_tibble(x))

df_pred_ibp <- bind_rows(list_kmn_pred_ibp) %>%
  add_column(N = L, t = 1:nrow(list_kmn_pred_ibp[[1]]),
             Model = "Gamma")

df_pred_ibp$N <- as.integer(df_pred_ibp$N)
df_pred_ibp <- df_pred_ibp %>%
  mutate( t = t + N )

joint_df_pred_bayes <- rbind(df_pred_poiss,df_pred_negbin,df_pred_ibp)

## 5.b) Good-Toulmin prediction
# compute
# labels_comb_ibp <- paste("N", Ns, sep = ".")
# 
# list_kmn_pred_gt <- vector(mode="list", length = length(Ns))
# names(list_kmn_pred_gt) <- labels_comb_ibp
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
#   list_kmn_pred_gt[[lab_comb]] <- predict_good_toulmin(N, M, sfs, cts, alternative = 0)$preds
#   
# }
# # Good-Toulmin predictions
# saveRDS(list_kmn_pred_gt, "mazziotta2016_application/plants/mazz_plants_op_gt_prediction.rds")

# read
list_kmn_pred_gt <- readRDS(file = "mazziotta2016_application/plants/mazz_plants_gt_prediction.rds")

list_kmn_pred_gt <- lapply(list_kmn_pred_gt, function(x) as_tibble(x))


df_pred_gt <- bind_rows(list_kmn_pred_gt) %>%
  add_column(N=L, t = 0:(nrow(list_kmn_pred_gt[[1]])-1),
             Model = "GT")

df_pred_gt$N <- as.integer(df_pred_gt$N)

df_pred_gt <- df_pred_gt %>%
  filter(t > N)


## 6.c) Chao estimates (for prediction)

list_chao_pred <- readRDS(file = "mazziotta2016_application/plants/mazz_plants_chao_rare.rds")
list_chao_pred <- lapply(list_chao_pred, function(x) as_tibble(x))

df_pred_chao <- bind_rows(list_chao_pred) %>%
  select(-c(lbs,ubs)) %>%
  rename(value = medians) %>%
  add_column(N=L, Model = "Chao")

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
joint_df_pred_bayes_plot <- joint_df_pred_bayes 
joint_df_pred_bayes_plot$Model <- factor(joint_df_pred_bayes_plot$Model, 
                                         levels = c("Gamma","Poisson","NegBin", "Chao", "GT"))

ggplot(joint_df_pred_bayes_plot, aes(x = t, y = medians, color = Model)) +
  geom_line(linetype = "dashed") +
  geom_ribbon(aes(ymin = lbs, ymax = ubs), alpha = 0.05) +
  geom_line(data = obs_sample, aes(t, obs), color="black", linetype="solid") +
  geom_line(data = df_pred_gt, aes(t, value)) +
  geom_line(data = df_pred_chao, aes(t, value)) +
  geom_vline(aes(xintercept = L), linetype="dashed", color = "grey") +
  facet_wrap(.~ N, scales = "free_x") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_tableau() +
  theme(aspect.ratio = 1)

ggsave(filename = "Plots_paper/plot_plants_op_prediction.png", width = 10, height = 4, dpi = 300, units = "in", device='png')








