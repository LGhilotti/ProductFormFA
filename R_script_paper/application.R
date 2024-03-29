rm(list=ls())

library(ProductFormFA)
library(xlsx)
library(readxl)
source("R_script_paper/Routine_Chao.R")
source("R_script_paper/utils.R")

# Choose dataset among the 3 available
type = "Fungi" # options: "Lichens", "Plants" 

data <- read_excel('R_script_paper/mazz2016_data.xls',
                   sheet = type) %>%
  filter(is.na(Species) == F) %>%
  column_to_rownames(var = "Species")

data <- t(data)
data <- data[, colSums(is.na(data))==0]
data <- data[, colSums(data)!=0]


# Number of sites and number of species
n <- nrow(data)
Kn <- ncol(data)
print(paste0("Number of sites: ", n ))
print(paste0("Number of species: ", Kn))

# Set extrapolation horizon 
M <- 1000

# Randomly reorder sites
seed <- 1234
set.seed(seed)

data_mat <- data[sample.int(n, size = n, replace = F),]

# Plot accumulation
plot(1:n, rarefaction(data_mat) )


# 2) Run the models on the data ----

# Empirical estimate of E(N) is obtained by Chiu
Nbar <- beta_binomial_estimator(data_mat)

# 2.1) Beta-Bernoulli with Poisson(lambda) mixture

# Initialization and MCMC setting 
init_PoissonBB <- list(alpha_0 = -1, s_0 = 1)
init_obj_PoissonBB <- initialization(model = "PoissonBB", init = init_PoissonBB )
mcmcparams_PoissonBB <- list(tau = 0.1, S = 3*10^4, n_burnin = 5*10^3, thin = 2)
mcmcparams_obj_PoissonBB <- mcmcparameters(model = "PoissonBB", mcmcparams = mcmcparams_PoissonBB)

# Hyperparameters elicitation 
hyper_PoissonBB <- list(a_alpha = 1, b_alpha = 0.1,
                        a_s = 2, b_s = 0.2,
                        lambda = Nbar)
prior_obj_PoissonBB <- prior(model = "PoissonBB", hyper = hyper_PoissonBB) 

# Fit the model
PoissonBB_fit <- GibbsFA(feature_matrix = data_mat, 
                         model = "PoissonBB", 
                         prior = prior_obj_PoissonBB, 
                         initialization = init_obj_PoissonBB, 
                         mcmcparams = mcmcparams_obj_PoissonBB)


# 2.2) Beta-Bernoulli with NB(n0, mu0) mixture

# Initialization and MCMC setting
init_NegBinBB <- list(alpha_0 = -1, s_0 = 1)
init_obj_NegBinBB <- initialization(model = "NegBinBB", init = init_NegBinBB )
mcmcparams_NegBinBB <- list(tau = 0.1, S = 3*10^4, n_burnin = 5*10^3, thin = 2)
mcmcparams_obj_NegBinBB <- mcmcparameters(model = "NegBinBB", mcmcparams = mcmcparams_NegBinBB)

# Hyperparameters elicitation 
c_fr <- 10

hyper_NegBinBB <- list(a_alpha = 1, b_alpha = 0.1,
                       a_s = 2, b_s = 0.2,
                       n0 = Nbar/(c_fr - 1), # n0, mu0 are set s.t. E(N) = Nbar, Var(N) = c_fr*E(N)
                       mu0 = 1/c_fr)
prior_obj_NegBinBB <- prior(model = "NegBinBB", hyper = hyper_NegBinBB) 


# Fit the model
NegBinBB_fit <- GibbsFA(feature_matrix = data_mat, 
                        model = "NegBinBB", 
                        prior = prior_obj_NegBinBB, 
                        initialization = init_obj_NegBinBB, 
                        mcmcparams = mcmcparams_obj_NegBinBB)


# 3) Estimate the total richness ----

PoissonBB_rich <- total_richness(object = PoissonBB_fit)
NegBinBB_rich <- total_richness(object = NegBinBB_fit)

# 4) Extrapolation ----

PoissonBB_extr <- extrapolation(object = PoissonBB_fit, M = M) 
NegBinBB_extr <- extrapolation(object = NegBinBB_fit, M = M) 

# 5) Rarefaction ----

PoissonBB_rare <- rarefaction(object = PoissonBB_fit) 
NegBinBB_rare <- rarefaction(object = NegBinBB_fit) 

# 6) Competitors: Chao and GT ------

# 6.1) Chao: competitor for rarefaction and extrapolation

# Determine the frequency vector of the training sets
Q_vec <- colSums(data_mat)
Q_vec <- Q_vec[Q_vec>0]

# Compute the curves with confidence intervals
fit_Chao <- iNEXT.Sam(Spec = Q_vec, T = n, endpoint = n + M)

rare_Chao <- as_tibble(fit_Chao[["q=0"]]) %>%
  select(-Cov.hat) %>%
  rename(medians = D0.hat, lbs = Norm.CI.Low, ubs = Norm.CI.High)

Chao_rare_extr <- as.data.frame(rare_Chao)


# 6.2) Smoothed Good-Toulmin: competitor for extrapolation

# Compute SFS vector and CTS vector
sfs <- tabulate(colSums(data_mat))
cts <- sapply(2:n, function(i) ncol(data_mat[1:i,colSums(data_mat[1:i,]) > 0])   )
cts <- c(0, sum(data_mat[1,]) , cts)

GT_extr <- predict_good_toulmin(n, M, sfs, cts, alternative = 0)$preds



# 7) Plot Richness: whole distributions for Poisson and NegBin

richness_PoissonBB_long <- as_tibble(PoissonBB_rich) %>%
  add_column(model= "Poisson")
      
richness_NegBinBB_long <- as_tibble(NegBinBB_rich) %>%
  add_column(model= "NegBin")

joint_richness_long <- bind_rows(richness_PoissonBB_long, richness_NegBinBB_long) %>%
  mutate(model = fct_relevel(model, c( "NegBin", "Poisson")))


ggplot(joint_richness_long, aes(x = value, color = model)) +
  stat_density(aes(x=value, colour=model),
               geom="line",position="identity", adjust = 2.5) +
  theme_light() +
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  xlab("# distinct features") + rremove("ylab") +
  scale_color_tableau() +
  theme(aspect.ratio = 1)



# 8) Plot Rarefaction: PoissonBB, NegBinBB and Chao

# PoissonBB
df_rare_PoissonBB <- as_tibble(t(bind_rows(as.data.frame(lapply(PoissonBB_rare, quantile, prob = c(0.025, 0.975))),
                               as.data.frame(lapply(PoissonBB_rare, mean))))) 
colnames(df_rare_PoissonBB) <- c("lbs", "ubs", "means")
df_rare_PoissonBB <- df_rare_PoissonBB %>%
  add_column(t = 1:nrow(df_rare_PoissonBB),
             model = "Poisson") %>%
  add_row(means = 0, ubs = 0, lbs = 0, t=0, model = "Poisson")


# NegBinBB
df_rare_NegBinBB <- as_tibble(t(bind_rows(as.data.frame(lapply(NegBinBB_rare, quantile, prob = c(0.025, 0.975))),
                                         as.data.frame(lapply(NegBinBB_rare, mean))))) 
colnames(df_rare_NegBinBB) <- c("lbs", "ubs", "means")
df_rare_NegBinBB <- df_rare_NegBinBB %>%
  add_column(t = 1:nrow(df_rare_NegBinBB),
             model = "NegBin") %>%
  add_row(means = 0, ubs = 0, lbs = 0, t=0, model = "NegBin")


# Chao
df_rare_Chao <- Chao_rare_extr %>%
  select(-c(lbs,ubs)) %>%
  rename(value = medians) %>%
  add_column(model = "Chao") %>%
  add_row(value = 0, t=0, model = "Chao") %>%
  filter(t <= n)


joint_df_rare_bayes <- rbind(df_rare_PoissonBB,df_rare_NegBinBB)

# plot
ggplot(joint_df_rare_bayes, aes(x = t, y = means, color = model)) +
  geom_line(linetype = "dashed", linewidth = 0.8) +
  geom_ribbon(aes(ymin = lbs, ymax = ubs), linewidth = 0.8, alpha = 0.1) +
  geom_line(data = df_rare_Chao, aes(t, value), linewidth = 0.8) +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_tableau() +
  theme(aspect.ratio = 1)




# 9) Plot Extrapolation: PoissonBB, NegBinBB, Chao and GT

# PoissonBB
df_extr_PoissonBB <- as_tibble(t(bind_rows(as.data.frame(lapply(PoissonBB_extr, quantile, prob = c(0.025, 0.975))),
                                         as.data.frame(lapply(PoissonBB_extr, mean))))) 
colnames(df_extr_PoissonBB) <- c("lbs", "ubs", "means")
df_extr_PoissonBB <- df_extr_PoissonBB %>%
  add_column(t = 1:nrow(df_extr_PoissonBB),
             model = "Poisson") %>%
  mutate( t = t + n ) %>%
  add_row(means = Kn, ubs = Kn, lbs = Kn, t=n, model = "Poisson")


# NegBinBB
df_extr_NegBinBB <- as_tibble(t(bind_rows(as.data.frame(lapply(NegBinBB_extr, quantile, prob = c(0.025, 0.975))),
                                        as.data.frame(lapply(NegBinBB_extr, mean))))) 
colnames(df_extr_NegBinBB) <- c("lbs", "ubs", "means")
df_extr_NegBinBB <- df_extr_NegBinBB %>%
  add_column(t = 1:nrow(df_extr_NegBinBB),
             model = "NegBin") %>%
  mutate( t = t + n ) %>%
  add_row(means = Kn, ubs = Kn, lbs = Kn, t=n, model = "NegBin")


# Chao
df_extr_Chao <- Chao_rare_extr %>%
  select(-c(lbs,ubs)) %>%
  rename(value = medians) %>%
  add_column(model = "Chao") %>%
  filter(t >= n)

# Good-Toulmin
df_extr_GT <- as_tibble(GT_extr) %>%
  add_column(t = 0:(length(GT_extr)-1),
             model = "GT") %>%
  filter(t >= n)


# Accumulation curve
accum <- rarefaction(data_mat)
accum_df <- data.frame("accum" = c(0, accum),
                       "t" = 0:length(accum))


joint_df_extr_bayes <- rbind(df_extr_PoissonBB,df_extr_NegBinBB)

# plot
ggplot(joint_df_extr_bayes, aes(x = t, y = means, color = model)) +
  geom_line(linetype = "dashed", linewidth = 0.8) +
  geom_ribbon(aes(ymin = lbs, ymax = ubs), linewidth = 0.8, alpha = 0.1) +
  geom_line(data = accum_df, aes(t, accum), color="black", linetype="solid") +
  geom_line(data = df_extr_Chao, aes(t, value), linewidth = 0.8) +
  geom_line(data = df_extr_GT, aes(t, value), linewidth = 0.8) +
  geom_vline(aes(xintercept = n), linetype="dashed", color = "grey") +
  xlab("# observations") + ylab("# distinct features") + 
  theme_light() + 
  theme(legend.position = "top") +
  scale_y_continuous(breaks = pretty_breaks()) +
  scale_x_continuous(breaks = pretty_breaks()) +
  scale_color_tableau() +
  theme(aspect.ratio = 1)








