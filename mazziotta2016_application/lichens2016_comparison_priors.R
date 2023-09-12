##############################################################
############ REAL DATA APPLICATION: lichens2016 ###################
#############################################################
rm(list=ls())

library(ProductFormFA)
library(readxl)
library(tidyverse)

lichens2016_mat <- read_excel('mazziotta2016_results/mazz2016_data.xls',
                              sheet = "Lichens") %>%
  filter(is.na(Species) == F) %>%
  column_to_rownames(var = "Species")

lichens2016_mat <- t(lichens2016_mat)

freq <- colSums(lichens2016_mat)
# to check that all species in the dataset are actually present 
sum(freq ==0) # must be =0

# number of sites
L = nrow(lichens2016_mat)

set.seed(12345)
lichens2016_mat <- lichens2016_mat[sample.int(L, size = L, replace = F),]

lichens2016_list <- create_features_list(lichens2016_mat)
plot_trajectory(lichens2016_list)

data_mat <- lichens2016_mat
data_list <- lichens2016_list

############ 1) Set parameters for the 3 models ###############

########### 1.1) Set parameters for BB (Poisson and Neg-Bin)

# Set prior hyperparameters specific for poisson
a_l_poiss <- 50
b_l_poiss <- 0.1
print(paste0("E(lambda) = ", a_l_poiss/ b_l_poiss))
print(paste0("Var(lambda) = ", a_l_poiss/ (b_l_poiss^2)))

# Set prior hyperparameters specific for NB
q_star_nb <- 0.01
print(paste0("E(nstar) = ", 1/ q_star_nb))
print(paste0("Var(nstar) = ", (1-q_star_nb)/ (q_star_nb^2)) )
alpha_p_nb <- 0.1
beta_p_nb <- 0.1
print(paste0("E(p) = ", alpha_p_nb/(alpha_p_nb + beta_p_nb)))
print(paste0("Var(p) = ",(alpha_p_nb*beta_p_nb)/((alpha_p_nb + beta_p_nb)^2 *(alpha_p_nb+beta_p_nb+1)) ))

# Set other hyperparameters
a_s_bb <- 2
b_s_bb <- 0.2
print(paste0("E(s) = ", a_s_bb/ b_s_bb))
print(paste0("Var(s) = ", a_s_bb/ (b_s_bb^2)))

# set 2 different possibilities for alpha 
a_alpha_bb_vague <- 0.7
b_alpha_bb_vague <- 0.07
print(paste0("E(alpha_bar) = ", a_alpha_bb_vague/ b_alpha_bb_vague))
print(paste0("Var(alpha_bar) = ", a_alpha_bb_vague/ (b_alpha_bb_vague^2)))

a_alpha_bb_conc <- 1
b_alpha_bb_conc <- 0.1
print(paste0("E(alpha_bar) = ", a_alpha_bb_conc/ b_alpha_bb_conc))
print(paste0("Var(alpha_bar) = ", a_alpha_bb_conc/ (b_alpha_bb_conc^2)))

# Set initial values for poisson 
lambda_0_poiss <- 100
# Set initial values for NB
nstar_0_nb <- 100
p_0_nb <- 0.2
# Set initial values for other parameters
alpha_bar_0_bb <- 1
s_0_bb <- 1


########## 1.2) Set parameters for the Gamma IBP and SB-SP

# Set prior hyperparameters for Gamma IBP
p_ibp <- 0.05
print(paste0("E(a) = ", 1/ p_ibp))
print(paste0("Var(a) = ", (1-p_ibp)/ (p_ibp^2)) )
r_ibp <- 1
t_ibp <- 0.1
print(paste0("E(b) = ", r_ibp/ t_ibp))
print(paste0("Var(b) = ", r_ibp/ (t_ibp^2)))
a_alpha_ibp <- 2
b_alpha_ibp <- 2
print(paste0("E(alpha) = ", a_alpha_ibp/ (a_alpha_ibp + b_alpha_ibp)))
print(paste0("Var(alpha) = ", a_alpha_ibp*b_alpha_ibp/ (a_alpha_ibp + b_alpha_ibp)^2 /(a_alpha_ibp+b_alpha_ibp+1)))
a_s_ibp <- 2
b_s_ibp <- 0.2
print(paste0("E(s) = ", a_s_ibp/ b_s_ibp))
print(paste0("Var(s) = ", a_s_ibp/ (b_s_ibp^2)))

print(paste0("E(theta) = ", a_s_ibp/ b_s_ibp -  a_alpha_ibp/ (a_alpha_ibp + b_alpha_ibp)))
print(paste0("Var(theta) = ", a_s_ibp/ (b_s_ibp^2) + 
               a_alpha_ibp*b_alpha_ibp/ (a_alpha_ibp + b_alpha_ibp)^2 /(a_alpha_ibp+b_alpha_ibp+1)))


# Set prior hyperparameters for SB-SP
p_sp <- 0.08
print(paste0("E(c) = ", (1-p_sp)/ p_sp))
print(paste0("Var(c) = ", (1-p_sp)/ (p_sp^2)) )
r_sp <- 0.1
t_sp <- 0.01
print(paste0("E(beta) = ", r_sp/ t_sp))
print(paste0("Var(beta) = ", r_sp/ (t_sp^2)))
a_alpha_sp <- 2
b_alpha_sp <- 2
print(paste0("E(alpha) = ", a_alpha_sp/ (a_alpha_sp + b_alpha_sp)))
print(paste0("Var(alpha) = ", a_alpha_sp*b_alpha_sp/ (a_alpha_sp + b_alpha_sp)^2 /(a_alpha_sp+b_alpha_sp+1)))


# Set initial values for the parameters of Gamma IBP
a_0_ibp <- 5
b_0_ibp <- 1
alpha_0_ibp <- 0.5
s_0_ibp <- 15

# Set initial values for the parameters of SB-SP
c_0_sp <- 10
beta_0_sp <- 10
alpha_0_sp <- 0.5


########## Set MCMC parameters (common to all 3 models): 
# SB-SP has same chain settings than Gamma IBP

S_poiss <- S_negbin <- S_ibp <-  8*10^4
n_burnin_poiss <- n_burnin_negbin <- n_burnin_ibp<-  10^4
thin_poiss <- thin_negbin <- thin_ibp <- 10
seed <- 12345
number_saved_iterations_poiss <- (S_poiss - n_burnin_poiss)/thin_poiss
number_saved_iterations_negbin <- (S_negbin - n_burnin_negbin)/thin_negbin
number_saved_iterations_ibp <- (S_ibp - n_burnin_ibp)/thin_ibp


###### 2) Run the algorithms ###############
gg_ntilde_poiss_vague <- data.frame(matrix(nrow = number_saved_iterations_poiss, ncol = 1))
colnames(gg_ntilde_poiss_vague) <- paste("N", L, sep = ".")
gg_ntilde_negbin_vague <- data.frame(matrix(nrow = number_saved_iterations_negbin, ncol = 1))
colnames(gg_ntilde_negbin_vague) <- paste("N", L, sep = ".")

gg_ntilde_poiss_conc <- data.frame(matrix(nrow = number_saved_iterations_poiss, ncol = 1))
colnames(gg_ntilde_poiss_conc) <- paste("N", L, sep = ".")
gg_ntilde_negbin_conc <- data.frame(matrix(nrow = number_saved_iterations_negbin, ncol = 1))
colnames(gg_ntilde_negbin_conc) <- paste("N", L, sep = ".")

params_poiss_vague <- data.frame(matrix(nrow = number_saved_iterations_poiss, ncol = 3))
colnames(params_poiss_vague) <- c("lambda", "alpha", "theta")
params_negbin_vague <- data.frame(matrix(nrow = number_saved_iterations_negbin, ncol = 4))
colnames(params_negbin_vague) <- c("nstar", "p", "alpha", "theta")

params_poiss_conc <- data.frame(matrix(nrow = number_saved_iterations_poiss, ncol = 3))
colnames(params_poiss_conc) <- c("lambda", "alpha", "theta")
params_negbin_conc <- data.frame(matrix(nrow = number_saved_iterations_negbin, ncol = 4))
colnames(params_negbin_conc) <- c("nstar", "p", "alpha", "theta")

################# 3) Run BB + Poisson - vague and conc ###########

#Set tau for MALA
tau_poiss <- 0.007

# vague
output_poiss_vague <- gibbs_sampler_poiss(Z = data_mat,
                                    lambda_0_poiss, alpha_bar_0_bb, s_0_bb,
                                    a_l_poiss, b_l_poiss, a_alpha_bb_vague, b_alpha_bb_vague, a_s_bb, b_s_bb,
                                    tau_poiss,
                                    S_poiss, n_burnin_poiss, thin_poiss, seed)

n_saved_iter_poiss_vague <- length(output_poiss_vague$lambda_vec)
lambda_chain_poiss_vague <- output_poiss_vague$lambda_vec
s_chain_poiss_vague <- output_poiss_vague$s_vec
alpha_bar_chain_poiss_vague <- output_poiss_vague$alpha_bar_vec
alpha_chain_poiss_vague <- - alpha_bar_chain_poiss_vague
theta_chain_poiss_vague <- s_chain_poiss_vague + alpha_bar_chain_poiss_vague

params_poiss_vague[["lambda"]] <- lambda_chain_poiss_vague
params_poiss_vague[["alpha"]] <- alpha_chain_poiss_vague
params_poiss_vague[["theta"]] <- theta_chain_poiss_vague

# conc
output_poiss_conc <- gibbs_sampler_poiss(Z = data_mat,
                                          lambda_0_poiss, alpha_bar_0_bb, s_0_bb,
                                          a_l_poiss, b_l_poiss, a_alpha_bb_conc, b_alpha_bb_conc, a_s_bb, b_s_bb,
                                          tau_poiss,
                                          S_poiss, n_burnin_poiss, thin_poiss, seed)

n_saved_iter_poiss_conc <- length(output_poiss_conc$lambda_vec)
lambda_chain_poiss_conc <- output_poiss_conc$lambda_vec
s_chain_poiss_conc <- output_poiss_conc$s_vec
alpha_bar_chain_poiss_conc <- output_poiss_conc$alpha_bar_vec
alpha_chain_poiss_conc <- - alpha_bar_chain_poiss_conc
theta_chain_poiss_conc <- s_chain_poiss_conc + alpha_bar_chain_poiss_conc

params_poiss_conc[["lambda"]] <- lambda_chain_poiss_conc
params_poiss_conc[["alpha"]] <- alpha_chain_poiss_conc
params_poiss_conc[["theta"]] <- theta_chain_poiss_conc

####### 4) Run BB + Negative-Binomial - vague and conc ################

# Set tau for MALA
tau_nb <- 0.01

# vague
output_negbin_vague <- gibbs_sampler_negbin_geometric(Z = data_mat,
                                                nstar_0_nb, p_0_nb,  s_0 = 5, alpha_bar_0 = 0.1,
                                                q_star_nb, alpha_p_nb, beta_p_nb,
                                                a_alpha_bb_vague, b_alpha_bb_vague, a_s_bb, b_s_bb,
                                                tau_nb, fixed = c(F,F,F,F),
                                                S_negbin, n_burnin_negbin, thin_negbin, seed)

n_saved_iter_negbin_vague <- length(output_negbin_vague$nstar_vec)
nstar_chain_negbin_vague <- output_negbin_vague$nstar_vec
p_chain_negbin_vague <- output_negbin_vague$p_vec
s_chain_negbin_vague <- output_negbin_vague$s_vec
alpha_bar_chain_negbin_vague <- output_negbin_vague$alpha_bar_vec
alpha_chain_negbin_vague <- - alpha_bar_chain_negbin_vague
theta_chain_negbin_vague <- s_chain_negbin_vague + alpha_bar_chain_negbin_vague

params_negbin_vague[["nstar"]] <- nstar_chain_negbin_vague
params_negbin_vague[["p"]] <- p_chain_negbin_vague
params_negbin_vague[["alpha"]] <- alpha_chain_negbin_vague
params_negbin_vague[["theta"]] <- theta_chain_negbin_vague

# Set tau for MALA
tau_nb <- 0.03

# conc
output_negbin_conc <- gibbs_sampler_negbin_geometric(Z = data_mat,
                                                      nstar_0_nb, p_0_nb,  s_0 = 5, alpha_bar_0 = 0.1,
                                                      q_star_nb, alpha_p_nb, beta_p_nb,
                                                      a_alpha_bb_conc, b_alpha_bb_conc, a_s_bb, b_s_bb,
                                                      tau_nb, fixed = c(F,F,F,F),
                                                      S_negbin, n_burnin_negbin, thin_negbin, seed)

n_saved_iter_negbin_conc <- length(output_negbin_conc$nstar_vec)
nstar_chain_negbin_conc <- output_negbin_conc$nstar_vec
p_chain_negbin_conc <- output_negbin_conc$p_vec
s_chain_negbin_conc <- output_negbin_conc$s_vec
alpha_bar_chain_negbin_conc <- output_negbin_conc$alpha_bar_vec
alpha_chain_negbin_conc <- - alpha_bar_chain_negbin_conc
theta_chain_negbin_conc <- s_chain_negbin_conc + alpha_bar_chain_negbin_conc

params_negbin_conc[["nstar"]] <- nstar_chain_negbin_conc
params_negbin_conc[["p"]] <- p_chain_negbin_conc
params_negbin_conc[["alpha"]] <- alpha_chain_negbin_conc
params_negbin_conc[["theta"]] <- theta_chain_negbin_conc

####### 6) Estimate limit distributions (Poiss/NB) ################
Kn = ncol(data_mat[,colSums(data_mat) > 0])

# Poisson
ntilde_chain_poiss_vague <- generate_Ntilde_chain_poiss(lambda_chain_poiss_vague,
                                                        alpha_chain_poiss_vague,
                                                  theta_chain_poiss_vague, n = L, Kn)
ntilde_chain_poiss_conc <- generate_Ntilde_chain_poiss(lambda_chain_poiss_conc,
                                                        alpha_chain_poiss_conc,
                                                        theta_chain_poiss_conc, n = L, Kn)

# Neg-Bin
ntilde_chain_negbin_vague <- generate_Ntilde_chain_negbin(nstar_chain_negbin_vague,
                                                          p_chain_negbin_vague,
                                                    alpha_chain_negbin_vague,
                                                    theta_chain_negbin_vague, n = L, Kn)

ntilde_chain_negbin_conc <- generate_Ntilde_chain_negbin(nstar_chain_negbin_conc,
                                                          p_chain_negbin_conc,
                                                          alpha_chain_negbin_conc,
                                                          theta_chain_negbin_conc, n = L, Kn)

gg_ntilde_poiss_vague[[paste0("N.",L)]] <- ntilde_chain_poiss_vague
gg_ntilde_poiss_conc[[paste0("N.",L)]] <- ntilde_chain_poiss_conc

gg_ntilde_negbin_vague[[paste0("N.",L)]] <- ntilde_chain_negbin_vague
gg_ntilde_negbin_conc[[paste0("N.",L)]] <- ntilde_chain_negbin_conc


############# 8) Save results:  MCMC convergence ####################

save(params_poiss_vague, file = "mazziotta2016_results/mazz_lichens_comp_params_poiss_vague.Rda")
save(params_poiss_conc, file = "mazziotta2016_results/mazz_lichens_comp_params_poiss_conc.Rda")
save(params_negbin_vague, file = "mazziotta2016_results/mazz_lichens_comp_params_negbin_vague.Rda")
save(params_negbin_conc, file = "mazziotta2016_results/mazz_lichens_comp_params_negbin_conc.Rda")

############ 9) Save results: samples from limiting distributions (Poiss/NB) ##############
# Poisson
save(gg_ntilde_poiss_vague, file = "mazziotta2016_results/mazz_lichens_comp_ntilde_poiss_vague.Rda")
save(gg_ntilde_poiss_conc, file = "mazziotta2016_results/mazz_lichens_comp_ntilde_poiss_conc.Rda")
# Negative Binomial
save(gg_ntilde_negbin_vague, file = "mazziotta2016_results/mazz_lichens_comp_ntilde_negbin_vague.Rda")
save(gg_ntilde_negbin_conc, file = "mazziotta2016_results/mazz_lichens_comp_ntilde_negbin_conc.Rda")




#####################################################################
####### ANALYSE RESULTS #############################################
####################################################################

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

###### 1) Read results:  MCMC convergence ####################
load(file = "mazziotta2016_results/mazz_lichens_comp_params_poiss_vague.Rda")
load(file = "mazziotta2016_results/mazz_lichens_comp_params_poiss_conc.Rda")
load(file =  "mazziotta2016_results/mazz_lichens_comp_params_negbin_vague.Rda")
load(file = "mazziotta2016_results/mazz_lichens_comp_params_negbin_conc.Rda")


lambda_chain_poiss_vague <- params_poiss_vague[["lambda"]]  
alpha_chain_poiss_vague <- params_poiss_vague[["alpha"]]
theta_chain_poiss_vague <- params_poiss_vague[["theta"]]

nstar_chain_negbin_vague <- params_negbin_vague[["nstar"]] 
p_chain_negbin_vague <- params_negbin_vague[["p"]]
alpha_chain_negbin_vague <- params_negbin_vague[["alpha"]]
theta_chain_negbin_vague <- params_negbin_vague[["theta"]]


lambda_chain_poiss_conc <- params_poiss_conc[["lambda"]]  
alpha_chain_poiss_conc <- params_poiss_conc[["alpha"]]
theta_chain_poiss_conc <- params_poiss_conc[["theta"]]

nstar_chain_negbin_conc <- params_negbin_conc[["nstar"]] 
p_chain_negbin_conc <- params_negbin_conc[["p"]]
alpha_chain_negbin_conc <- params_negbin_conc[["alpha"]]
theta_chain_negbin_conc <- params_negbin_conc[["theta"]]

###### 2) Read results: samples from limiting distributions (Poiss/NB) ##############
load(file = "mazziotta2016_results/mazz_lichens_comp_ntilde_poiss_vague.Rda")
load(file = "mazziotta2016_results/mazz_lichens_comp_ntilde_poiss_conc.Rda")
load(file = "mazziotta2016_results/mazz_lichens_comp_ntilde_negbin_vague.Rda")
load(file = "mazziotta2016_results/mazz_lichens_comp_ntilde_negbin_conc.Rda")

###### 4) Read the data ###############################
data_mat <- readRDS(file = "mazziotta2016_results/mazz_lichens_data_mat.rds")
data_list <- create_features_list(data_mat)
L <- nrow(data_mat)


###### 5) Check MCMC convergence###################

###### check mcmc mixing poisson - vague
samples_poiss_vague <- mcmc.list(mcmc(params_poiss_vague))
samples_ggs_poiss_vague <- ggs(samples_poiss_vague, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_poiss_vague) 
ggs_running(samples_ggs_poiss_vague)

###### check mcmc mixing poisson - conc
samples_poiss_conc <- mcmc.list(mcmc(params_poiss_conc))
samples_ggs_poiss_conc <- ggs(samples_poiss_conc, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_poiss_conc) 
ggs_running(samples_ggs_poiss_conc)

###### check mcmc mixing negbin - vague
samples_negbin_vague <- mcmc.list(mcmc(params_negbin_vague))
samples_ggs_negbin_vague <- ggs(samples_negbin_vague, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_negbin_vague) 
ggs_running(samples_ggs_negbin_vague)

###### check mcmc mixing negbin - conc
samples_negbin_conc <- mcmc.list(mcmc(params_negbin_conc))
samples_ggs_negbin_conc <- ggs(samples_negbin_conc, keep_original_order = TRUE)
ggs_traceplot(samples_ggs_negbin_conc) 
ggs_running(samples_ggs_negbin_conc)

###### 8) Richness - vague
gg_ntilde_poiss_vague_long <- gg_ntilde_poiss_vague %>%
  rename(estimate = paste0("N.", L)) %>%
  add_column(Model= "BBmixP-vague") %>%
  filter(estimate < 3000)
gg_ntilde_poiss_conc_long <- gg_ntilde_poiss_conc %>%
  rename(estimate = paste0("N.", L)) %>%
  add_column(Model= "BBmixP-conc")%>%
  filter(estimate < 3000)

gg_ntilde_negbin_vague_long <- gg_ntilde_negbin_vague %>%
  rename(estimate = paste0("N.", L)) %>%
  add_column(Model= "BBmixNB-vague") %>%
  filter(estimate < 3000)
gg_ntilde_negbin_conc_long <- gg_ntilde_negbin_conc %>%
  rename(estimate = paste0("N.", L)) %>%
  add_column(Model= "BBmixNB-conc") %>%
  filter(estimate < 3000)

joint_total_long <- bind_rows(gg_ntilde_poiss_vague_long, 
                              gg_ntilde_poiss_conc_long,
                              gg_ntilde_negbin_vague_long,
                              gg_ntilde_negbin_conc_long) %>%
  mutate(Model = fct_relevel(Model, c("BBmixP-vague","BBmixP-conc",
                                      "BBmixNB-vague", "BBmixNB-conc"))) 

# plot
ggplot(joint_total_long, aes(x = estimate, color = Model)) + 
  #geom_density(alpha = 0.5, linewidth = 0.8) +
  stat_density(aes(x=estimate, colour=Model),
               geom="line",position="identity") +
  geom_vline(xintercept = ncol(data_mat), color="grey", linetype="dashed", linewidth=0.8) + # observed features in the whole dataset
  #facet_wrap(~N, labeller = labeller(N = ~ paste("n = ", .x)), scales = "free", nrow = 1) +
  theme_bw() + 
  theme(strip.text.x = element_text(size = 12))+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  xlab("# distinct features") + rremove("ylab") #+
  # scale_color_manual(values = c("BBmixP" = "forestgreen",
  #                               "BBmixNB" = "royalblue1")) 

ggsave(filename = "Plots_paper/plot_lichens2016_comp_richness.png", width = 8, height = 5, dpi = 300, units = "in", device='png')
