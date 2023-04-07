#########################################################
############ Script for some checks ###################
#######################################################

#install.packages("comprehenr")
library(comprehenr)

###############################################################
################ BB with Negative-Binomial(n*,p) #############
##############################################################
set.seed(1234)

# Credible intervals for the number of features (not conditional)
ci_kmn_negbin_bb <- CI_Kmn_negbin_BB(alpha = -1, theta = 10, m = 100, n = 0,
                                     Kn = 0, nstar = 50, p = 0.5, lev = 0.95)

# mean of the number of features on the whole sample
ci_kmn_negbin_bb$means[100]

# Generate from buffet procedure from beginning
# (returns the features and the number of new features for new customers)
buff_negbin_bb <- buffet_negbin_BB(alpha = -1, theta = 10, n = 10000, nstar = 100, p = 0.8)
n_feat <- length(buff_negbin_bb$counts)

# Empirical Bayes via EFPF maximization (n*,p parametrization)
eb_est_fixed <- EB_EFPF_fixed_negbin_BB(n = length(buff_negbin_bb$num_new), 
                                        counts = buff_negbin_bb$counts, 
                                        pars_0 = c(-10, 15, 50, 0.2),
                                        opt = c(T,T,T,T))

# objective function as function of alpha
alpha_grid = seq(from = -3, to = -0.01, length.out = 1000)
theta = 10; nstar= 100; p=0.8; alpha = -1;

ev_nlEFPF_alpha <- to_vec(for(alp in alpha_grid) 
  neg_log_EFPF_negbin_BB_rep_all(pars = c(alp,alp+theta,nstar,p),
                                 n = length(buff_negbin_bb$num_new),
                                 counts = buff_negbin_bb$counts))

plot(alpha_grid, ev_nlEFPF_alpha, pch=16, cex=0.5)

# objective function as function of theta
theta_grid = seq(from = 2, to = 50, length.out = 1000)
alpha = -1; nstar= 100; p=0.8; 

ev_nlEFPF_theta <- to_vec(for(theta in theta_grid) 
  neg_log_EFPF_negbin_BB_rep_all(pars = c(alpha,alpha+theta,nstar,p),
                                 n = length(buff_negbin_bb$num_new),
                                 counts = buff_negbin_bb$counts))

plot(theta_grid, ev_nlEFPF_theta, pch=16, cex=0.5)

# objective function as function of nstar
nstar_grid = seq(from = 10, to = 100, length.out = 1000)
alpha=-1; theta = 10; p=0.8; 

ev_nlEFPF_nstar <- to_vec(for(nstar in nstar_grid) 
  neg_log_EFPF_negbin_BB_rep_all(pars = c(alpha,alpha+theta,nstar,p),
                                 n = length(buff_negbin_bb$num_new),
                                 counts = buff_negbin_bb$counts))

plot(nstar_grid, ev_nlEFPF_nstar, pch=16, cex=0.5)

# objective function as function of p
p_grid = seq(from = 0.01, to = 0.99, length.out = 1000)
alpha=-1; theta = 10; nstar=100; 

ev_nlEFPF_p <- to_vec(for(p in p_grid) 
  neg_log_EFPF_negbin_BB_rep_all(pars = c(alpha,alpha+theta,nstar,p),
                                 n = length(buff_negbin_bb$num_new),
                                 counts = buff_negbin_bb$counts))

plot(p_grid, ev_nlEFPF_p, pch=16, cex=0.5)

##if n large <- alpha_min si avvicina a 0 e 
##if n small <- marginalmente hanno tutte minimo nel punto corretto, perchè non ritrova?

########################################
###### plot jointly nstar and p ########
########################################
#install.packages("plotly")
library(plotly)

# vectorize the "neg_log_EFPF_negbin_BB_rep_all" to take nstar_grid and p_grid 
# as vector
vec_nlEFPF_nstar_p <- Vectorize(function(nstar, p, alpha, theta, n, counts){
  return (neg_log_EFPF_negbin_BB_rep_all(pars = c(alpha,alpha+theta,nstar,p),
                                         n = n,
                                         counts = counts))
}, vectorize.args=c("nstar", "p"))

eval_nlEFPF_nstar_p <- outer(nstar_grid,p_grid, vec_nlEFPF_nstar_p, alpha=alpha, theta=theta, 
      n= length(buff_negbin_bb$num_new), counts = buff_negbin_bb$counts)

fig <- plot_ly(x=~nstar_grid, y=~p_grid, z=~eval_nlEFPF_nstar_p)%>% add_surface(
  contours = list(
    z = list(
      show=TRUE,
      usecolormap=TRUE,
      highlightcolor="#ff0000",
      project=list(z=TRUE)
    )
  )
)
fig <- fig %>% layout(
  scene = list(
    camera=list(
      eye = list(x=1.87, y=0.88, z=-0.64)
    )
  )
)

fig


##########################################################
###### Mean-Variance Negative-Binomial parametrization ###
##########################################################
# Empirical Bayes via EFPF maximization (m-v parametrization)
eb_est_mv <- EB_EFPF_fixed_negbin_mv_BB(n = length(buff_negbin_bb$num_new), 
                                        counts = buff_negbin_bb$counts, 
                                        pars_0 = c(-10, 15, 50, 100), # (alpha,theta,mu,sigma2)
                                        opt = c(T,T,T,T))

eb_est_np <- eb_est_mv
eb_est_np[3] <- (eb_est_mv[3]**2)/(eb_est_mv[4] - eb_est_mv[3])
eb_est_np[4] <- eb_est_mv[3]/eb_est_mv[4]

# objective function as function of alpha
mu <- nstar*(1-p)/p
sigma2 <- nstar*(1-p)/(p**2)

ev_nlEFPF_alpha_mv <- to_vec(for(alp in alpha_grid) 
  neg_log_EFPF_negbin_mv_BB_rep_all(pars = c(alp,alp+theta,mu,sigma2),
                                 n = length(buff_negbin_bb$num_new),
                                 counts = buff_negbin_bb$counts))

plot(alpha_grid, ev_nlEFPF_alpha_mv, pch=16, cex=0.5)

# objective function as function of theta
mu <- nstar*(1-p)/p
sigma2 <- nstar*(1-p)/(p**2)

ev_nlEFPF_theta_mv <- to_vec(for(th in theta_grid) 
  neg_log_EFPF_negbin_mv_BB_rep_all(pars = c(alpha,alpha+th,mu,sigma2),
                                    n = length(buff_negbin_bb$num_new),
                                    counts = buff_negbin_bb$counts))

plot(theta_grid, ev_nlEFPF_theta_mv, pch=16, cex=0.5)


########################################
###### plot jointly nstar and p ########
########################################

mu_grid <- seq(from = 0.1, to = 100, length.out = 100)
t_grid <- seq(from = 0.1, to = 100, length.out = 100)

# vectorize the "neg_log_EFPF_negbin_BB_rep_all" to take nstar_grid and p_grid 
# as vector
vec_nlEFPF_mu_t <- Vectorize(function(mu, t, alpha, theta, n, counts){
  return (neg_log_EFPF_negbin_mv_BB_rep_all(pars = c(alpha,alpha+theta,mu,t),
                                         n = n,
                                         counts = counts))
}, vectorize.args=c("mu", "t"))

eval_nlEFPF_mu_t <- outer(mu_grid,t_grid, vec_nlEFPF_mu_t, alpha=alpha, theta=theta, 
                             n= length(buff_negbin_bb$num_new), counts = buff_negbin_bb$counts)

fig2 <- plot_ly(x=~mu_grid, y=~t_grid, z=~eval_nlEFPF_mu_t)%>% add_surface(
  contours = list(
    z = list(
      show=TRUE,
      usecolormap=TRUE,
      highlightcolor="#ff0000",
      project=list(z=TRUE)
    )
  )
)
fig2 <- fig2 %>% layout(
  scene = list(
    camera=list(
      eye = list(x=1.87, y=0.88, z=-0.64)
    )
  )
)

fig2

#####################################################################
##### try to fix the variance hyperparameter of the NB ###############
######################################################################
true_mean <- nstar*(1-p)/p
true_var <- nstar*(1-p)/(p**2)
# Empirical Bayes via EFPF maximization (m-v parametrization)
eb_est_fixed_mv <- EB_EFPF_fixed_negbin_mv_BB(n = length(buff_negbin_bb$num_new), 
                                              counts = buff_negbin_bb$counts, 
                                              pars_0 = c(-1, 15, 10, true_var),
                                              opt = c(T,T,T,F))

eb_est_fixed_np <- eb_est_fixed_mv
eb_est_fixed_np[3] <- (eb_est_fixed_mv[3]**2)/(eb_est_fixed_mv[4] - eb_est_fixed_mv[3])
eb_est_fixed_np[4] <- eb_est_fixed_mv[3]/eb_est_fixed_mv[4]
eb_est_fixed_np
# theta alto fa si che le possibili feature N siano tutte usate, cosi N può 
# essere centrato intorno a k. Oppure anche mettere N


#################################################################################
### CHECK IF THE MMLE RECOVERS WHEN REPLICATES OF THE POPULATION IS AVAILABLE ###
#################################################################################

set.seed(1234)
n_replicates <- 10
l_counts_negbin_bb <- vector("list",n_replicates) # list of counts of features for each replicate
# alpha = -100; theta = 101; n = 10000; nstar = 100; p = 0.8
alpha = -100; theta = 101; n = 100000; nstar = 200; p = 0.8

true_mean <- nstar*(1-p)/p
true_mean

true_var <- nstar*(1-p)/(p**2)
true_var

# Generate from buffet procedure multiple times and store in list
for (i in 1:n_replicates){
  
  buff_negbin_bb <- buffet_negbin_BB(alpha, theta, n, nstar, p)
  l_counts_negbin_bb[[i]] <- buff_negbin_bb$counts
}


# -> with replicates
avg_n_obs_feat <- length(unlist(l_counts_negbin_bb))/n_replicates
avg_n_obs_feat
eb_efpf_mv_replicates <- EB_EFPF_fixed_negbin_mv_BB_replicates(n, 
                                                               counts = l_counts_negbin_bb,
                                                               pars_0 = c(-1, 5, 10, 20),
                                                               opt = c(T,T,T,T))
eb_efpf_mv_replicates
# -> single realization of the sample
eb_efpf_mv_single <- EB_EFPF_fixed_negbin_mv_BB(n,
                                                counts = l_counts_negbin_bb[[1]],
                                                pars_0 = c(-3, 10, 30, 60),
                                                opt = c(T,T,T,T))
eb_efpf_mv_single



#############################################################
#### CHECK ON THE BUFFET PROCESS ############################
##########################################################
printList <- function(list, max) {
  
  for (item in 1:max) {
    
    print(list[[item]])
    
  }
}


alpha = -100; theta = 101; n = 10000; nstar = 100; p = 0.8


set.seed(1234)
buff_negbin_bb_slow <- buffet_negbin_BB(alpha, theta, n, nstar, p)

