rm(list=ls())
library(ProductFormFA)
library(tidyverse)
library(scales)
library(ggpubr)
library(ggthemes)
library(patchwork)
library(latex2exp)
library(dplyr, warn.conflicts = FALSE)

# Load the Work space
mechanism <- "custom"
load(paste0("R_script_tommaso/eb_EFPF_",mechanism,"_fit_estimate_singledataset.RData"))
Kn <- sapply(Ns, function(n) sum(colSums(data_mat[1:n,]) > 0)  )


#### Richness -------

# 1) Plot Richness: Point plot expected value (number of features)
labels_comb_bb <- paste(rep(paste("n_train", Ns, sep = "."), each = length(Nbars)+1),
                        c("Nbar.emp" , paste("Nbar", Nbars, sep = ".")), sep=":")

# PoissonBB
richness_EFPF_PoissonBB_df <- tibble(estimate = unname(sapply(list_eb_EFPF_fit_PoissonBB, function(x)
  total_richness(x)$lambda_post + ncol(x$feature_matrix)) ) ) %>%
  add_column(Model = "Poisson BB", 
             Nbar = rep(c("EB", Nbars), length(Ns)),
             n_train = rep(Ns, each = length(Nbars)+ 1),
             n_train_idx = rep(c(1,2,3), each = length(Nbars)+ 1)) 

# NegBinBB
richness_EFPF_NegBinBB_df <- tibble(estimate = numeric(), Model = character(),
                                    Nbar = character(), 
                                    n_train = integer(), n_train_idx = integer())

for (var_fct_NegBinBB in vars_fct_NegBinBB){
  list_eb_EFPF_fit_NegBinBB_var <- list_eb_EFPF_fit_NegBinBB[[paste0("var_fct.", var_fct_NegBinBB)]]
  
  richness_EFPF_NegBinBB_df_var <- tibble(estimate = unname(sapply(list_eb_EFPF_fit_NegBinBB_var, function(x)
    total_richness(x)$mu0_post + ncol(x$feature_matrix)) ) ) %>%
    add_column(Model = paste0("NegBinomial BB x",var_fct_NegBinBB) , 
               Nbar = rep(c("EB", Nbars), length(Ns)),
               n_train = rep(Ns, each = length(Nbars)+ 1),
               n_train_idx = rep(c(1,2,3), each = length(Nbars)+ 1))
  
  richness_EFPF_NegBinBB_df <- bind_rows(richness_EFPF_NegBinBB_df, richness_EFPF_NegBinBB_df_var)
  
}


joint_richness_long <- bind_rows(richness_EFPF_PoissonBB_df,
                                 richness_EFPF_NegBinBB_df) %>%
  mutate(Model = fct_relevel(Model, c("Poisson BB", 
                                      paste0("NegBinomial BB x", vars_fct_NegBinBB))) )

# n_train.labs <- paste0("n = ", Ns,", Kn = ", Kn )
# names(n_train.labs) <- Ns

joint_richness_long <- joint_richness_long %>%
  mutate(n_train_latex = paste0(r"($n = $)", n_train, r"($, K_n =$)", Kn[n_train_idx]))


# plots estimator
ggplot(joint_richness_long, aes( y=estimate, x=Model, shape = Nbar)) +
  geom_point(size = 2) +
  facet_wrap(~ TeX(n_train_latex, output = "character"),
             labeller = label_parsed, #custom_labeller,
             scales = "free_x", nrow = 1) +
  theme_light() +
  geom_hline(aes(yintercept = H), linetype = "dashed") +
  theme(legend.position = "top") +
  ylab("Posterior mean of N") +
  scale_y_continuous(breaks = pretty_breaks()) +
  rremove("xlab") +
  scale_shape_discrete(name = "Prior mean of N") +
  theme(aspect.ratio = 1) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

