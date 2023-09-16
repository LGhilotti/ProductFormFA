####
##### MODEL 3 ############################
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


###### 2) Read results: samples from limiting distributions (Poiss/NB) ##############
load(file = "chiu_model_simulation/m3/m3_ntilde_poiss_large.Rda")
load(file = "chiu_model_simulation/m3/m3_ntilde_negbin_large.Rda")
#load(file = "chiu_model_simulation/m3/m3_ntilde_negbin_prior.Rda")



###### 4) Read the data ###############################
data_mat <- readRDS(file = "chiu_model_simulation/m3/m3_data_mat_large.rds")
data_list <- create_features_list(data_mat)
L <- nrow(data_mat)
Ns <- c(800,2000, 10000)

Nbars <- c(200)
c_fr <- 10



####### 6) Plot limiting distributions (Poiss/NB) ##############
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

ggsave(filename = "Plots_paper/plot_m3_richness.png", width = 10, height = 4, dpi = 300, units = "in", device='png')




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
  
  # # Chao2 estimator
  # table_richness <- table_richness %>%
  #   add_row(N = N, Nbar = NA, Model = "Chao2",
  #           estimator = chao2_estimator(train_mat),
  #           Model_spec = "Chao2")
  # 
  # # Improved Chao2 estimator
  # table_richness <- table_richness %>%
  #   add_row(N = N, Nbar = NA, Model = "ImpChao2",
  #           estimator = improved_chao2_estimator(train_mat),
  #           Model_spec = "ImpChao2")
  # 
  
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



