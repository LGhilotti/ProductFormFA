##############################################################
############ REAL DATA APPLICATION: plants2016 ###################
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

# ##########################################################
# #### (A) Training/test approach ############################
# ##########################################################
# 
# ###### 1) Read results:  MCMC convergence ####################
# load(file = "mazziotta2016_results/mazz_plants_params_poiss.Rda")
# load(file =  "mazziotta2016_results/mazz_plants_params_negbin.Rda")
# load(file =  "mazziotta2016_results/mazz_plants_params_ibp.Rda" )
# load(file =  "mazziotta2016_results/mazz_plants_params_sp.Rda" )
# 
# ###### 2) Read results: samples from limiting distributions (Poiss/NB) ##############
# load(file = "mazziotta2016_results/mazz_plants_ntilde_poiss.Rda")
# load(file = "mazziotta2016_results/mazz_plants_ntilde_negbin.Rda")
# 
# ###### 3) Read results: CI for extrapolation (Poiss/NB/Gamma) ################
# list_kmn_pred_test_poiss <- readRDS(file = "mazziotta2016_results/mazz_plants_ci_poiss.rds")
# list_kmn_pred_test_negbin <- readRDS(file = "mazziotta2016_results/mazz_plants_ci_negbin.rds")
# list_kmn_pred_test_ibp <- readRDS(file = "mazziotta2016_results/mazz_plants_ci_ibp.rds")
# list_kmn_pred_test_sp <- readRDS(file = "mazziotta2016_results/mazz_plants_ci_sp.rds")
# 
# ###### 3.b) Read results: CI for insample (Poiss/NB/Gamma) ################
# list_kn_rarefaction_poiss <- readRDS(file = "mazziotta2016_results/mazz_plants_ci_insample_poiss.rds")
# list_kn_rarefaction_negbin <- readRDS(file = "mazziotta2016_results/mazz_plants_ci_insample_negbin.rds")
# list_kn_rarefaction_ibp <- readRDS(file = "mazziotta2016_results/mazz_plants_ci_insample_ibp.rds")
# list_kn_rarefaction_sp <- readRDS(file = "mazziotta2016_results/mazz_plants_ci_insample_sp.rds")
# 
# 
# ###### 4) Read the data ###############################
# data_mat <- readRDS(file = "mazziotta2016_results/mazz_plants_data_mat.rds")
# data_list <- create_features_list(data_mat)
# L <- nrow(data_mat)
# Ms <- sapply(list_kmn_pred_test_poiss, function(l) length(l$medians))
# Ns <- L - Ms
# 
# 
# ###### 5) Check MCMC convergence###################
# 
# ###### check mcmc mixing poisson
# samples_poiss <- mcmc.list(mcmc(params_poiss))
# samples_ggs_poiss <- ggs(samples_poiss, keep_original_order = TRUE)
# ggs_traceplot(samples_ggs_poiss) + 
#   facet_wrap(~Parameter, nrow = length(Ns), scales = "free")
# 
# ###### check mcmc mixing negbin
# samples_negbin <- mcmc.list(mcmc(params_negbin))
# samples_ggs_negbin <- ggs(samples_negbin, keep_original_order = TRUE)
# ggs_traceplot(samples_ggs_negbin) + 
#   facet_wrap(~Parameter, nrow = length(Ns), scales = "free")
# 
# ###### check mcmc mixing Gamma ibp
# samples_ibp <- mcmc.list(mcmc(params_ibp))
# samples_ggs_ibp <- ggs(samples_ibp, keep_original_order = TRUE)
# ggs_traceplot(samples_ggs_ibp) + 
#   facet_wrap(~Parameter, nrow = length(Ns), scales = "free")
# 
# ###### check mcmc mixing SB-SP
# samples_sp <- mcmc.list(mcmc(params_sp))
# samples_ggs_sp <- ggs(samples_sp, keep_original_order = TRUE)
# ggs_traceplot(samples_ggs_sp) + 
#   facet_wrap(~Parameter, nrow = length(Ns), scales = "free")
# 
# 
# ###### 6) Rarefaction curves
# gg_kn_rarefaction_all  <- vector(mode="list", length = length(Ns))
# 
# for (j in 1:length(Ns)){
#   N <- Ns[j]
#   
#   train_mat <- data_mat[1:N,]
#   # convert the binary matrix into list of features
#   train_list <- create_features_list(train_mat)
#   
#   gg_kn_rarefaction_all[[j]] <- plot_Kn_median_and_rarefaction_all(
#     train_list = train_list,
#     ci_poiss = list_kn_rarefaction_poiss[[paste0("N.",N)]], 
#     ci_negbin = list_kn_rarefaction_negbin[[paste0("N.",N)]],
#     ci_ibp = list_kn_rarefaction_ibp[[paste0("N.",N)]],
#     n_avg = 100) +
#     ggtitle(paste0("n = ", N)) + theme(plot.title = element_text(size=12)) 
#   
#   
#   if (j != 1){
#     gg_kn_rarefaction_all[[j]] <- gg_kn_rarefaction_all[[j]] +
#       theme(#axis.text.y = element_blank(),
#         #axis.ticks.y = element_blank(),
#         axis.title.y = element_blank() )
#   }
#   
# }
# 
# # Print plots
# fig_plants2016_rare_ <- wrap_plots(gg_kn_rarefaction_all, nrow = 1, ncol = 3) +   plot_layout(guides = "collect") & theme(legend.position = 'right') & xlab(NULL) & theme(plot.margin = margin(5.5, 5.5, 0, 5.5))
# 
# fig_plants2016_rare <- wrap_elements(panel = fig_plants2016_rare_) +
#   labs(tag = "# observations") +
#   theme(
#     plot.tag = element_text(size = rel(1)),
#     plot.tag.position = "bottom") 
# 
# ggsave(filename = "Plots_paper/plot_plants2016_rare.png", width = 10, height = 4, dpi = 300, units = "in", device='png')
# 
# 
# ###### 7) Prediction
# gg_kmn_pred_test_all  <- vector(mode="list", length = length(Ns))
# 
# for (j in 1:length(Ns)){
#   N <- Ns[j]
#   M <- L - N
#   
#   train_mat <- data_mat[1:N,]
#   test_mat <- data_mat[(N+1):L, ]
#   # convert the binary matrix into list of features
#   train_list <- create_features_list(train_mat)
#   test_list <- create_features_list(test_mat)
#   
#   
#   gg_kmn_pred_test_all[[j]] <- plot_Kmn_median_pred_and_test_all(
#     train_list = train_list,
#     test_list = test_list,
#     ci_poiss = list_kmn_pred_test_poiss[[paste0("N.",N)]], 
#     ci_negbin = list_kmn_pred_test_negbin[[paste0("N.",N)]],
#     ci_ibp = list_kmn_pred_test_ibp[[paste0("N.",N)]],
#     n_avg = 100) +
#     ggtitle(paste0("n = ", N)) + theme(plot.title = element_text(size=12))
#   
#   
#   if (j != 1){
#     gg_kmn_pred_test_all[[j]] <- gg_kmn_pred_test_all[[j]] +
#       theme(#axis.text.y = element_blank(),
#         #axis.ticks.y = element_blank(),
#         axis.title.y = element_blank() )
#   }
#   
# }
# 
# # Print plots
# fig_plants2016_pred_ <- wrap_plots(gg_kmn_pred_test_all, nrow = 1, ncol = 3) +   plot_layout(guides = "collect") & theme(legend.position = 'right') & xlab(NULL) & theme(plot.margin = margin(5.5, 5.5, 0, 5.5))
# 
# fig_plants2016_pred <- wrap_elements(panel = fig_plants2016_pred_) +
#   labs(tag = "# observations") +
#   theme(
#     plot.tag = element_text(size = rel(1)),
#     plot.tag.position = "bottom") 
# 
# ggsave(filename = "Plots_paper/plot_plants2016_pred.png", width = 10, height = 4, dpi = 300, units = "in", device='png')
# 
# 
# ###### 8) Richness
# gg_ntilde_poiss_long <- gather(gg_ntilde_poiss, training, estimate, 
#                                paste0("N.",Ns[1]):paste0("N.",Ns[length(Ns)]), 
#                                factor_key=TRUE) %>%
#   add_column(Model= "BBmixP", N = rep(Ns, each = nrow(gg_ntilde_poiss))) 
# 
# gg_ntilde_negbin_long <- gather(gg_ntilde_negbin, training, estimate, 
#                                 paste0("N.",Ns[1]):paste0("N.",Ns[length(Ns)]), 
#                                 factor_key=TRUE) %>%
#   add_column(Model= "BBmixNB", N = rep(Ns, each = nrow(gg_ntilde_negbin)))
# 
# 
# joint_total_long <- bind_rows(gg_ntilde_poiss_long, gg_ntilde_negbin_long) %>%
#   mutate(Model = fct_relevel(Model, c("BBmixP", "BBmixNB"))) 
# 
# # plot
# ggplot(joint_total_long, aes(x = estimate, color = Model)) + 
#   #geom_density(alpha = 0.5, linewidth = 0.8) +
#   stat_density(aes(x=estimate, colour=Model), bw=10,
#                geom="line",position="identity") +
#   geom_vline(xintercept = ncol(data_mat), color="grey", linetype="dashed", linewidth=0.8) + # observed features in the whole dataset
#   facet_wrap(~N, labeller = labeller(N = ~ paste("n = ", .x)), scales = "free", nrow = 1) +
#   theme_bw() + 
#   theme(strip.text.x = element_text(size = 12))+
#   theme(#panel.grid.major = element_blank(),
#     #panel.grid.minor = element_blank(),
#     strip.background = element_blank(),
#     panel.border = element_rect(colour = "black", fill=NA)) +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   scale_y_continuous(breaks = pretty_breaks()) +
#   #xlim(c(100,180)) +
#   xlab("# distinct features") + rremove("ylab") +
#   scale_color_manual(values = c("BBmixP" = "forestgreen",
#                                 "BBmixNB" = "royalblue1")) 
# 
# ggsave(filename = "Plots_paper/plot_plants2016_richness.png", width = 10, height = 4, dpi = 300, units = "in", device='png')



##########################################################
#### (B) Only prediction #################################
##########################################################
rm(list=ls())
library(ggmcmc)
library(coda)
library(ggpubr)
library(scales)
library(tidyverse)
library(gridExtra)
library(grid)
library(patchwork)

##### Single dataset -> Ntilde (Poiss/NB) and extrapolation (Poiss/NB/Gamma) #####

###### 1) Read results:  MCMC convergence ####################
load(file = "mazziotta2016_application/mazz_plants_op_params_poiss.Rda")
load(file =  "mazziotta2016_application/mazz_plants_op_params_negbin.Rda")
#load(file =  "mazziotta2016_application/mazz_plants_op_params_negbin_prior.Rda")
load(file =  "mazziotta2016_application/mazz_plants_op_params_ibp.Rda" )
#load(file =  "mazziotta2016_application/mazz_plants_op_params_sp.Rda" )

###### 2) Read results: samples from limiting distributions (Poiss/NB) ##############
load(file = "mazziotta2016_application/mazz_plants_op_ntilde_poiss.Rda")
load(file = "mazziotta2016_application/mazz_plants_op_ntilde_negbin.Rda")
#load(file = "mazziotta2016_application/mazz_plants_op_ntilde_negbin_prior.Rda")

###### 3) Read results: CI for extrapolation (Poiss/NB/Gamma) ################
list_kmn_pred_test_poiss <- readRDS(file = "mazziotta2016_application/mazz_plants_op_ci_pred_poiss.rds")
list_kmn_pred_test_negbin <- readRDS(file = "mazziotta2016_application/mazz_plants_op_ci_pred_negbin.rds")
#list_kmn_pred_test_negbin_prior <- readRDS(file = "mazziotta2016_application/mazz_plants_op_ci_pred_negbin_prior.rds")
list_kmn_pred_test_ibp <- readRDS(file = "mazziotta2016_application/mazz_plants_op_ci_pred_ibp.rds")
#list_kmn_pred_test_sp <- readRDS(file = "mazziotta2016_application/mazz_plants_op_ci_sp.rds")

# ###### 3.b) Read results: CI for insample (Poiss/NB/Gamma) ################
# list_kn_rarefaction_poiss <- readRDS(file = "mazziotta2016_application/mazz_plants_op_ci_insample_poiss.rds")
# list_kn_rarefaction_negbin <- readRDS(file = "mazziotta2016_application/mazz_plants_op_ci_insample_negbin.rds")
# #list_kn_rarefaction_negbin_prior <- readRDS(file = "mazziotta2016_application/mazz_plants_op_ci_insample_negbin_prior.rds")
# list_kn_rarefaction_ibp <- readRDS(file = "mazziotta2016_application/mazz_plants_op_ci_insample_ibp.rds")
# #list_kn_rarefaction_sp <- readRDS(file = "mazziotta2016_application/mazz_plants_op_ci_insample_sp.rds")


###### 4) Read the data ###############################
data_mat <- readRDS(file = "mazziotta2016_application/mazz_plants_data_mat.rds")
data_list <- create_features_list(data_mat)
L <- nrow(data_mat)

Nbars <- c(200,400,600)

###### 5) Check MCMC convergence###################

# Poisson 

###### check mcmc mixing poisson Nbar = 200
params_poiss_ <- params_poiss[c(paste0("alpha:Nbar.", Nbars[1] ), 
                                paste0("theta:Nbar.", Nbars[1])) ]
samples_poiss <- mcmc.list(mcmc(params_poiss_))
samples_ggs_poiss <- ggs(samples_poiss, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_poiss) + 
  facet_wrap(~Parameter, nrow = length(Nbars), scales = "free")

effectiveSize(params_poiss_)

###### check mcmc mixing poisson Nbar = 400
params_poiss_ <- params_poiss[c(paste0("alpha:Nbar.", Nbars[2]), 
                                paste0("theta:Nbar.", Nbars[2])) ]
samples_poiss <- mcmc.list(mcmc(params_poiss_))
samples_ggs_poiss <- ggs(samples_poiss, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_poiss) + 
  facet_wrap(~Parameter, nrow = length(Nbars), scales = "free")

effectiveSize(params_poiss_)

###### check mcmc mixing poisson Nbar = 600
params_poiss_ <- params_poiss[c(paste0("alpha:Nbar.", Nbars[3]), 
                                paste0("theta:Nbar.",  Nbars[3])) ]
samples_poiss <- mcmc.list(mcmc(params_poiss_))
samples_ggs_poiss <- ggs(samples_poiss, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_poiss) + 
  facet_wrap(~Parameter, nrow = length(Nbars), scales = "free")

effectiveSize(params_poiss_)

# Negative Binomial

###### check mcmc mixing negbin (fixed) Nbar = 200
params_negbin_ <- params_negbin[c(paste0("alpha:Nbar.", Nbars[1]), 
                                  paste0("theta:Nbar.", Nbars[1])) ]
samples_negbin <- mcmc.list(mcmc(params_negbin_))
samples_ggs_negbin <- ggs(samples_negbin, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_negbin) + 
  facet_wrap(~Parameter, nrow = length(Nbars), scales = "free")

effectiveSize(params_negbin_)

###### check mcmc mixing negbin (fixed) Nbar = 400
params_negbin_ <- params_negbin[c(paste0("alpha:Nbar.", Nbars[2]), 
                                  paste0("theta:Nbar.", Nbars[2])) ]
samples_negbin <- mcmc.list(mcmc(params_negbin_))
samples_ggs_negbin <- ggs(samples_negbin, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_negbin) + 
  facet_wrap(~Parameter, nrow = length(Nbars), scales = "free")

effectiveSize(params_negbin_)

###### check mcmc mixing negbin (fixed) Nbar = 600
params_negbin_ <- params_negbin[c(paste0("alpha:Nbar.", Nbars[3]), 
                                  paste0("theta:Nbar.", Nbars[3])) ]
samples_negbin <- mcmc.list(mcmc(params_negbin_))
samples_ggs_negbin <- ggs(samples_negbin, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_negbin) + 
  facet_wrap(~Parameter, nrow = length(Nbars), scales = "free")

effectiveSize(params_negbin_)

# Gamma ibp

###### check mcmc mixing Gamma ibp Nbar=200
params_ibp_ <- params_ibp[c(paste0("alpha:Nbar.", Nbars[1]), 
                            paste0("theta:Nbar.", Nbars[1])) ]
samples_ibp <- mcmc.list(mcmc(params_ibp_))
samples_ggs_ibp <- ggs(samples_ibp, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_ibp) + 
  facet_wrap(~Parameter, nrow = length(Nbars), scales = "free")

effectiveSize(params_ibp_)

###### check mcmc mixing Gamma ibp Nbar=400
params_ibp_ <- params_ibp[c(paste0("alpha:Nbar.", Nbars[2]), 
                            paste0("theta:Nbar.", Nbars[2])) ]
samples_ibp <- mcmc.list(mcmc(params_ibp_))
samples_ggs_ibp <- ggs(samples_ibp, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_ibp) + 
  facet_wrap(~Parameter, nrow = length(Nbars), scales = "free")

effectiveSize(params_ibp_)

###### check mcmc mixing Gamma ibp Nbar=600
params_ibp_ <- params_ibp[c(paste0("alpha:Nbar.", Nbars[3]), 
                            paste0("theta:Nbar.", Nbars[3])) ]
samples_ibp <- mcmc.list(mcmc(params_ibp_))
samples_ggs_ibp <- ggs(samples_ibp, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_ibp) + 
  facet_wrap(~Parameter, nrow = length(Nbars), scales = "free")

effectiveSize(params_ibp_)


####### 6) Plot limiting distributions (Poiss/NB) ##############
labels_comb <- paste("Nbar", Nbars, sep = ".")


gg_ntilde_poiss_long <- gather(gg_ntilde_poiss, training, estimate, 
                               labels_comb, 
                               factor_key=TRUE) %>%
  add_column(Model= "BBmixP", 
             Nbar = rep(Nbars, each = nrow(gg_ntilde_poiss)) )


gg_ntilde_negbin_long <- gather(gg_ntilde_negbin, training, estimate, 
                                labels_comb, 
                                factor_key=TRUE) %>%
  add_column(Model= "BBmixNB", 
             Nbar = rep(Nbars, each = nrow(gg_ntilde_negbin)) )


joint_total_long <- bind_rows(gg_ntilde_poiss_long, gg_ntilde_negbin_long) %>%
  mutate(Model = fct_relevel(Model, c("BBmixP", "BBmixNB"))) 

# plot
ggplot(joint_total_long, aes(x = estimate, color = Model)) + 
  #geom_density(alpha = 0.5, linewidth = 0.8) +
  stat_density(aes(x=estimate, colour=Model),
               geom="line",position="identity") +
  facet_wrap(~ Nbar, labeller = labeller(Nbar = ~ paste("Nbar = ", .x)), scales = "free") +
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

ggsave(filename = "Plots_paper/plot_plants2016_op_richness.png", width = 10, height = 4, dpi = 300, units = "in", device='png')


######## 7) Plot Extrapolation curve (Poiss/NB/Gamma) ################

gg_kmn_pred_test_all  <- vector(mode="list", length = length(Nbars))
names(gg_kmn_pred_test_all) <- labels_comb


for (v in 1:length(Nbars)){
  Nbar <- Nbars[v]
  
  lab_comb <- paste0("Nbar.",Nbar)
  
  gg_kmn_pred_test_all[[lab_comb]] <- plot_Kmn_median_pred_all(
    data_list = data_list,
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



# Print plots
fig_plants2016_op_pred_ <- wrap_plots(gg_kmn_pred_test_all, nrow = 1, ncol = 3) +   plot_layout(guides = "collect") & theme(legend.position = 'right') & xlab(NULL) & theme(plot.margin = margin(5.5, 5.5, 0, 5.5))

fig_plants2016_op_pred <- wrap_elements(panel = fig_plants2016_op_pred_) +
  labs(tag = "# observations") +
  theme(
    plot.tag = element_text(size = rel(1)),
    plot.tag.position = "bottom") 

ggsave(filename = "Plots_paper/plot_plants2016_op_pred.png", width = 10, height = 4, dpi = 300, units = "in", device='png')

