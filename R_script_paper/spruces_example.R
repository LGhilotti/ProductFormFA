#
# Application to spruces data ####
#

rm(list=ls())
library(ProductFormFA)
library(tidyverse)
library(scales)
library(ggpubr)
library(ggthemes)

source("R_script_paper/utils.R")

# Read the CSV file
data <- read.csv("R_script_paper/spruce_data_a_1_b_20.csv", header = FALSE)

data <- data[, colSums(is.na(data))==0]
data <- data[, colSums(data)!=0]

# Training set
n_train <- 60
data_mat <- as.matrix(data)[1:n_train, ]
data_mat <- data_mat[, colSums(data_mat)!=0]

# Number of sites and number of species
n <- nrow(data_mat) # coincides with n_train
Kn <- ncol(data_mat)
print(paste0("Number of sites: ", n ))
print(paste0("Number of species: ", Kn))


# Plot accumulation
accum_df <- tibble( x = 0:n,
                    n_feat = c(0,rarefaction(data_mat[1:n,], n_reorderings = 1)))

ggplot(accum_df, aes(x = x, y = n_feat)) +
  geom_point(color="black", shape = 19, size = 0.1) + 
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  theme(aspect.ratio = 1) + 
  scale_color_tableau()


# 1) Complete estimation ----------

# Initial parameters for optimization
eb_init_BB <- list(alpha = -10, s = 10, Nhat_prime = 300)
eb_known_BB <- list()

eb_params_obj_BB <- eb_params(model = "BB", 
                              init = eb_init_BB, known = eb_known_BB )

# PoissonBB
eb_EFPF_fit_PoissonBB <- GibbsFA_eb(feature_matrix = data_mat, 
                                    model = "PoissonBB", 
                                    type = "EFPF",
                                    eb_params =  eb_params_obj_BB)


###### 2) Oracle on the mark parameters ----------

# Initial parameters for optimization
eb_init_BB_oracle <- list(Nhat_prime = 300)
eb_known_BB_oracle <- list(alpha = -1, s = 20)

eb_params_obj_BB_oracle <- eb_params(model = "BB", 
                              init = eb_init_BB_oracle, known = eb_known_BB_oracle )

# PoissonBB
eb_EFPF_fit_PoissonBB_oracle <- GibbsFA_eb(feature_matrix = data_mat, 
                                    model = "PoissonBB", 
                                    type = "EFPF",
                                    eb_params =  eb_params_obj_BB_oracle)


