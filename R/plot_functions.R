#' Plot function for the credible intervals of Kmn
#'
#' This function allows to plot the credible intervals of Kmn
#'
#' @param ci [list] it contains means, upper-bounds and lower-bounds of the credible intervals
#'
#' @export
#' @import ggplot2
#' @import scales
#'
plot_Kmn <- function(ci){
  
  m <- length(ci$means)
  
  bands <- data.frame(
    x = 1:m,
    means = ci$means,
    lbs = ci$lbs,
    ubs = ci$ubs
  )
  
  ggplot(bands, aes(x, means)) + # ggplot2 plot with confidence intervals
    geom_line(col = "darkblue") +
    geom_ribbon(aes(ymin = lbs, ymax = ubs), alpha = 0.1, fill = "darkblue") +
    xlab("m") + ylab(expression(K[m]^n)) + theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_x_continuous(breaks = pretty_breaks())
}

#####################################################################
#' Plot function for the credible intervals of Kmn
#'
#' This function allows to plot the credible intervals of Kmn
#'
#' @param ci [list] it contains medians, upper-bounds and lower-bounds of the credible intervals
#'
#' @export
#' @import ggplot2
#' @import scales
#'
plot_Kmn_median <- function(ci){
  
  m <- length(ci$medians)
  
  bands <- data.frame(
    x = 1:m,
    medians = ci$medians,
    lbs = ci$lbs,
    ubs = ci$ubs
  )
  
  ggplot(bands, aes(x, medians)) + # ggplot2 plot with confidence intervals
    geom_line(col = "darkblue") +
    geom_ribbon(aes(ymin = lbs, ymax = ubs), alpha = 0.1, fill = "darkblue") +
    xlab("m") + ylab(expression(K[m]^n)) + theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_x_continuous(breaks = pretty_breaks())
}


##############################################################
#####################################################################
#' Plot function for the credible intervals of Kn
#'
#' This function allows to plot the credible intervals of Kn
#'
#' @param ci [list] it contains medians, upper-bounds and lower-bounds of the credible intervals
#'
#' @export
#' @import ggplot2
#' @import scales
#'
plot_Kn_median <- function(ci){
  
  m <- length(ci$medians)
  
  bands <- data.frame(
    x = 1:m,
    medians = ci$medians,
    lbs = ci$lbs,
    ubs = ci$ubs
  )
  
  ggplot(bands, aes(x, medians)) + # ggplot2 plot with confidence intervals
    geom_line(col = "darkblue") +
    geom_ribbon(aes(ymin = lbs, ymax = ubs), alpha = 0.1, fill = "darkblue") +
    xlab("m") + ylab(expression(K[n])) + theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_x_continuous(breaks = pretty_breaks())
}


##############################################################
#####################################################
#' Plot function for the credible intervals of Kmn, given initial sample
#'
#' This function allows to plot the credible intervals of Kmn given initial sample
#'
#' @param data_list [list] list of features in the initial sample 
#' @param ci [list] it contains means, upper-bounds and lower-bounds of the credible intervals
#'
#' @export
#' @import ggplot2
#' @import scales
#'
plot_Kn_median_and_sample <- function(data_list, ci){
  
  N <- length(ci$medians)
  
  cum_nfeat <- sapply(1:N, function(n) length(unique(unlist(data_list[1:n]))))
  
  obs_sample <- data.frame(
    x = 0:N,
    nfeat = c(0,cum_nfeat)
  )
  

  bands <- data.frame(
    x = 0:N,
    medians = c(0, ci$medians),
    lbs = c(0, ci$lbs),
    ubs = c(0, ci$ubs)
  )
  
  
  ggfig <- ggplot(bands, aes(x,medians) ) +
    geom_line(col = "red", linetype = "dashed") +
    geom_ribbon(aes(ymin = lbs, ymax = ubs), alpha = 0.1, fill = "red") +
    geom_line(data = obs_sample, aes(x, nfeat), color="black", linetype="solid", size=0.5) +
    xlab("n = # observations") + ylab(expression(K[n])) + theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("Kn within sample") +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_x_continuous(breaks = pretty_breaks())
  
  return (ggfig)
}

##############################################################
#####################################################
#' Plot function for the credible intervals of Kmn, given initial sample
#'
#' This function allows to plot the credible intervals of Kmn given initial sample
#'
#' @param data_list [list] list of features in the initial sample 
#' @param ci [list] it contains means, upper-bounds and lower-bounds of the credible intervals
#'
#' @export
#' @import ggplot2
#' @import scales
#'
plot_Kn_median_and_rarefaction <- function(train_list, ci, n_avg){
  
  N <- length(ci$medians)
  
  cum_nfeat <- rep(0, N)
  
  for (j in 1:n_avg){
    ord <- sample(1:N)
    cum_nfeat <- cum_nfeat + sapply(1:N, function(n) length(unique(unlist(train_list[ord][1:n]))))
  }
  
  cum_nfeat <- cum_nfeat/n_avg
  
  obs_sample <- data.frame(
    x = 0:N,
    nfeat = c(0,cum_nfeat)
  )
  
  
  bands <- data.frame(
    x = 0:N,
    medians = c(0, ci$medians),
    lbs = c(0, ci$lbs),
    ubs = c(0, ci$ubs)
  )
  
  
  ggfig <- ggplot(bands, aes(x,medians) ) +
    geom_line(col = "red", linetype = "dashed") +
    geom_ribbon(aes(ymin = lbs, ymax = ubs), alpha = 0.1, fill = "red") +
    geom_line(data = obs_sample, aes(x, nfeat), color="black", linetype="solid", size=0.5) +
    xlab("n = # observations") + ylab(expression(K[n])) + theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("Kn within sample and rarefaction curve") +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_x_continuous(breaks = pretty_breaks())
  
  return (ggfig)
}



#####################################################
#' Plot function for the credible intervals IN-SAMPLE of Kn, and the observed, 
#' for all the 3 models 
#' 
#'
#' This function allows to plot the credible intervals of Kn, a posteriori,
#' and the observed data as well
#'
#' @param train_list [list] 
#' @param ci_poiss [list] it contains means, upper-bounds and lower-bounds of the CI for Poisson
#' @param ci_negbin [list] it contains means, upper-bounds and lower-bounds of the CI for NB
#' @param ci_ibp [list] it contains means, upper-bounds and lower-bounds of the CI for IBP
#' @param n_avg [integer] 
#' 
#' @export
#' @import ggplot2
#' @import scales
#'
plot_Kn_median_and_rarefaction_all <- function(
    train_list, 
    ci_poiss, 
    ci_negbin,
    ci_ibp,
    n_avg){
  
  N <- length(ci_poiss$medians)
  
  cum_nfeat <- rep(0, N)
  
  for (j in 1:n_avg){
    ord <- sample(1:N)
    cum_nfeat <- cum_nfeat + sapply(1:N, function(n) length(unique(unlist(train_list[ord][1:n]))))
  }
  
  cum_nfeat <- cum_nfeat/n_avg
  
  obs_sample <- data.frame(
    x = 0:N,
    nfeat = c(0,cum_nfeat)
  )
  
  # bands
  bands_poiss <- data.frame(
    x = 0:N,
    medians = c(0, ci_poiss$medians),
    lbs = c(0, ci_poiss$lbs),
    ubs = c(0, ci_poiss$ubs),
    model = "Poiss"
  )
  
  bands_negbin <- data.frame(
    x = 0:N,
    medians = c(0, ci_negbin$medians),
    lbs = c(0, ci_negbin$lbs),
    ubs = c(0, ci_negbin$ubs),
    model = "NegBin"
  )
  
  bands_ibp <- data.frame(
    x = 0:N,
    medians = c(0, ci_ibp$medians),
    lbs = c(0, ci_ibp$lbs),
    ubs = c(0, ci_ibp$ubs),
    model = "IBP"
  )
  
  bands <- rbind(bands_poiss, bands_negbin, bands_ibp)
  
  ggfig <- ggplot(bands, aes(x,medians, color = model) ) +
    geom_line(linetype = "dashed") +
    geom_ribbon(aes(ymin = lbs, ymax = ubs, fill = model), alpha = 0.1) +
    geom_line(data = obs_sample, aes(x, nfeat), color="black", linetype="solid", linewidth=0.5) +
    xlab("# observations") + ylab("# distinct features") + theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_x_continuous(breaks = pretty_breaks())
  
  return (ggfig)
}


#####################################################
#' Plot function for the credible intervals of Kmn, given initial sample
#'
#' This function allows to plot the credible intervals of Kmn given initial sample
#'
#' @param N [integer] dimension of initial sample
#' @param train_list [list] list of features in the initial sample 
#' @param ci [list] it contains means, upper-bounds and lower-bounds of the credible intervals
#'
#' @export
#' @import ggplot2
#' @import scales
#'
plot_Kmn_given_sample <- function(N, train_list, ci){
  
  m <- length(ci$means)
  
  cum_nfeat <- sapply(1:N, function(n) length(unique(unlist(train_list[1:n]))))
    
  init_sample <- data.frame(
    x = 0:N,
    nfeat = c(0,cum_nfeat)
  )
  
  nfeat_sample <- cum_nfeat[N]
  
  bands <- data.frame(
    x = N:(N+m),
    means = c(nfeat_sample, ci$means + nfeat_sample),
    lbs = c(nfeat_sample, ci$lbs + nfeat_sample),
    ubs = c(nfeat_sample, ci$ubs + nfeat_sample)
  )
  
  
  ggfig <- ggplot(bands, aes(x,means) ) +
    geom_line(col = "red", linetype = "dashed") +
    geom_ribbon(aes(ymin = lbs, ymax = ubs), alpha = 0.1, fill = "red") +
    geom_line(data = init_sample, aes(x, nfeat), color="red", linetype="dashed") +
    geom_segment(aes(x = N, y = 0, xend = N, yend = ubs[m+1] + 10), color="grey",
                 linetype="dashed", size=1) +
    xlab("# observations") + ylab("# distinct features") + theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0("N = ", N)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_x_continuous(breaks = pretty_breaks())
  
  return (ggfig)
}

#####################################################
#' Plot function for the credible intervals of Kmn, given initial sample
#'
#' This function allows to plot the credible intervals of Kmn given initial sample
#'
#' @param N [integer] dimension of initial sample
#' @param train_list [list] list of features in the initial sample 
#' @param ci [list] it contains means, upper-bounds and lower-bounds of the credible intervals
#'
#' @export
#' @import ggplot2
#' @import scales
#'
plot_Kmn_median_pred_and_rarefaction <- function(train_list, ci, n_avg){
  
  m <- length(ci$medians)
  N <- length(train_list)
  
  cum_nfeat <- rep(0, N)
  
  for (j in 1:n_avg){
    ord <- sample(1:N)
    cum_nfeat <- cum_nfeat + sapply(1:N, function(n) length(unique(unlist(train_list[ord][1:n]))))
  }
  
  cum_nfeat <- cum_nfeat/n_avg
  
  init_sample <- data.frame(
    x = 0:N,
    nfeat = c(0,cum_nfeat)
  )
  
  nfeat_sample <- cum_nfeat[N]
  
  bands <- data.frame(
    x = N:(N+m),
    medians = c(nfeat_sample, ci$medians + nfeat_sample),
    lbs = c(nfeat_sample, ci$lbs + nfeat_sample),
    ubs = c(nfeat_sample, ci$ubs + nfeat_sample)
  )
  
  
  ggfig <- ggplot(bands, aes(x,medians) ) +
    geom_line(col = "red", linetype = "dashed") +
    geom_ribbon(aes(ymin = lbs, ymax = ubs), alpha = 0.1, fill = "red") +
    geom_line(data = init_sample, aes(x, nfeat), color="black", linetype="solid", size=0.5) +
    geom_segment(aes(x = N, y = 0, xend = N, yend = ubs[m+1] + 10), color="grey",
                 linetype="dashed", size=1) +
    xlab("# observations") + ylab("# distinct features") + theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0("N = ", N)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_x_continuous(breaks = pretty_breaks())
  
  return (ggfig)
}

#####################################################
#' Plot function for the credible intervals of Kmn, given initial sample, and the observed 
#' 
#'
#' This function allows to plot the credible intervals of Kmn given initial sample,
#' and the observed test as well
#'
#' @param N [integer] dimension of initial sample
#' @param data_list [list] list of features in the whole sample
#' @param ci [list] it contains means, upper-bounds and lower-bounds of the credible intervals
#'
#' @export
#' @import ggplot2
#' @import scales
#'
plot_Kmn_given_sample_with_observed <- function(N, data_list, ci){
  
  m <- length(ci$means)
  
  cum_nfeat <- sapply(1:(N+m), function(n) length(unique(unlist(data_list[1:n]))))
  
  obs_sample <- data.frame(
    x = 0:(N+m),
    nfeat = c(0,cum_nfeat)
  )
  
  nfeat_sample <- cum_nfeat[N]
  
  bands <- data.frame(
    x = N:(N+m),
    means = c(nfeat_sample, ci$means + nfeat_sample),
    lbs = c(nfeat_sample, ci$lbs + nfeat_sample),
    ubs = c(nfeat_sample, ci$ubs + nfeat_sample)
  )
  
  
  ggfig <- ggplot(bands, aes(x,means) ) +
    geom_line(col = "red", linetype = "dashed") +
    geom_ribbon(aes(ymin = lbs, ymax = ubs), alpha = 0.1, fill = "red") +
    geom_line(data = obs_sample, aes(x, nfeat), color="black", linetype="solid", size=0.5) +
    geom_segment(aes(x = N, y = 0, xend = N, yend = max(ubs[m+1],cum_nfeat[N+m]) + 10), color="grey",
                 linetype="dashed", size=1) +
    xlab("# observations") + ylab("# distinct features") + theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0("N = ", N)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_x_continuous(breaks = pretty_breaks())
  
  return (ggfig)
}

#####################################################
#' Plot function for the credible intervals of Kmn, given initial sample, and the observed 
#' 
#'
#' This function allows to plot the credible intervals of Kmn given initial sample,
#' and the observed test as well
#'
#' @param N [integer] dimension of initial sample
#' @param data_list [list] list of features in the whole sample
#' @param ci [list] it contains means, upper-bounds and lower-bounds of the credible intervals
#'
#' @export
#' @import ggplot2
#' @import scales
#'
plot_Kmn_median_pred_and_test <- function(train_list, test_list, ci, n_avg){
  
  m <- length(test_list)
  N <- length(train_list)
  L <- m + N
  
  # train set
  cum_nfeat <- rep(0, N)
  
  for (j in 1:n_avg){
    ord <- sample(1:N)
    cum_nfeat <- cum_nfeat + sapply(1:N, function(n) length(unique(unlist(train_list[ord][1:n]))))
  }
  
  cum_nfeat <- cum_nfeat/n_avg
  
  nfeat_sample <- cum_nfeat[N]
  
  names_train_features <- unique(unlist(train_list))
    
  # test set
  cum_nfeat_test <- rep(0, m)
  
  for (j in 1:n_avg){
    ord <- sample(1:m)
    cum_nfeat_test <- cum_nfeat_test + 
      sapply(1:m, function(n) length(setdiff(unique(unlist(test_list[ord][1:n])), 
                                     names_train_features))) 
  }
  
  cum_nfeat_test <- cum_nfeat_test/n_avg
  
  obs_sample <- data.frame(
    x = 0:L,
    nfeat = c(0, cum_nfeat, nfeat_sample + cum_nfeat_test)
  )
  
  # credible bands
  bands <- data.frame(
    x = N:(N+m),
    medians = c(nfeat_sample, ci$medians + nfeat_sample),
    lbs = c(nfeat_sample, ci$lbs + nfeat_sample),
    ubs = c(nfeat_sample, ci$ubs + nfeat_sample)
  )
  
  
  ggfig <- ggplot(bands, aes(x,medians) ) +
    geom_line(col = "red", linetype = "dashed") +
    geom_ribbon(aes(ymin = lbs, ymax = ubs), alpha = 0.1, fill = "red") +
    geom_line(data = obs_sample, aes(x, nfeat), color="black", linetype="solid", size=0.5) +
    geom_segment(aes(x = N, y = 0, xend = N, yend = max(ubs[m+1],nfeat_sample + cum_nfeat_test[m]) + 10), 
                 color="grey", linetype="dashed", size=1) +
    xlab("# observations") + ylab("# distinct features") + theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0("N = ", N)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_x_continuous(breaks = pretty_breaks())
  
  return (ggfig)
}


#####################################################
#' Plot function for the credible intervals of Kmn, given initial sample, and the observed, 
#' for all the 3 models 
#' 
#'
#' This function allows to plot the credible intervals of Kmn given initial sample,
#' and the observed test as well
#'
#' @param train_list [list] 
#' @param test_list [list] 
#' @param ci_poiss [list] it contains means, upper-bounds and lower-bounds of the CI for Poisson
#' @param ci_negbin [list] it contains means, upper-bounds and lower-bounds of the CI for NB
#' @param ci_ibp [list] it contains means, upper-bounds and lower-bounds of the CI for IBP
#' @param n_avg [integer] 
#' 
#' @export
#' @import ggplot2
#' @import scales
#'
plot_Kmn_median_pred_and_test_all <- function(
  train_list,
  test_list,
  ci_poiss, 
  ci_negbin,
  ci_ibp,
  n_avg){
  
  m <- length(test_list)
  N <- length(train_list)
  L <- m + N
  
  # train set
  cum_nfeat <- rep(0, N)
  
  for (j in 1:n_avg){
    ord <- sample(1:N)
    cum_nfeat <- cum_nfeat + sapply(1:N, function(n) length(unique(unlist(train_list[ord][1:n]))))
  }
  
  cum_nfeat <- cum_nfeat/n_avg
  
  nfeat_sample <- cum_nfeat[N]
  
  names_train_features <- unique(unlist(train_list))
  
  # test set
  cum_nfeat_test <- rep(0, m)
  
  for (j in 1:n_avg){
    ord <- sample(1:m)
    cum_nfeat_test <- cum_nfeat_test + 
      sapply(1:m, function(n) length(setdiff(unique(unlist(test_list[ord][1:n])), 
                                             names_train_features))) 
  }
  
  cum_nfeat_test <- cum_nfeat_test/n_avg
  
  obs_sample <- data.frame(
    x = 0:L,
    nfeat = c(0, cum_nfeat, nfeat_sample + cum_nfeat_test)
  )
  
  # credible bands Poiss
  bands_poiss <- data.frame(
    x = N:(N+m),
    medians = c(nfeat_sample, ci_poiss$medians + nfeat_sample),
    lbs = c(nfeat_sample, ci_poiss$lbs + nfeat_sample),
    ubs = c(nfeat_sample, ci_poiss$ubs + nfeat_sample), 
    model = "Poiss"
  )
  
  # credible bands NegBin
  bands_negbin <- data.frame(
    x = N:(N+m),
    medians = c(nfeat_sample, ci_negbin$medians + nfeat_sample),
    lbs = c(nfeat_sample, ci_negbin$lbs + nfeat_sample),
    ubs = c(nfeat_sample, ci_negbin$ubs + nfeat_sample), 
    model = "NegBin"
  )
  
  # credible bands IBP
  bands_ibp <- data.frame(
    x = N:(N+m),
    medians = c(nfeat_sample, ci_ibp$medians + nfeat_sample),
    lbs = c(nfeat_sample, ci_ibp$lbs + nfeat_sample),
    ubs = c(nfeat_sample, ci_ibp$ubs + nfeat_sample),
    model = "IBP"
  )
  
  bands <- rbind(bands_poiss, bands_negbin, bands_ibp)
    
  ggfig <- ggplot(bands, aes(x,medians, color = model) ) +
    geom_line(linetype = "dashed") +
    geom_ribbon(aes(ymin = lbs, ymax = ubs, fill = model), alpha = 0.1) +
    geom_line(data = obs_sample, aes(x, nfeat), color="black", linetype="solid", linewidth=0.5)  +
    geom_segment(aes(x = N, y = 0, xend = N, yend = max(max(ubs),nfeat_sample + cum_nfeat_test[m]) + 10), 
                 color="grey", linetype="dashed", size=1) +
    xlab("# observations") + ylab("# distinct features") + theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0("Extrapolated rarefaction curve, N = ", N)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_x_continuous(breaks = pretty_breaks())
  
  return (ggfig)
}



#####################################################
#' Plot function for the credible intervals of Kmn, given initial sample, for 
#' both efpf and MM, and the observed trajectory
#' 
#'
#' This function allows to plot the credible intervals of Kmn given initial sample,
#' for both efpf and MM, and the observed test as well
#'
#' @param N [integer] dimension of initial sample
#' @param data_list [list] list of features in the whole sample
#' @param ci_efpf [list] it contains means, upper-bounds and lower-bounds of the CI for efpf
#' @param ci_mm [list] it contains means, upper-bounds and lower-bounds of the CI for mm
#' 
#' @export
#' @import ggplot2
#' @import scales
#'
plot_Kmn_given_sample_with_observed_both <- function(N, data_list, ci_efpf, ci_mm){
  
  m <- length(ci_efpf$means)
  
  cum_nfeat <- sapply(1:(N+m), function(n) length(unique(unlist(data_list[1:n]))))
  
  obs_sample <- data.frame(
    x = 0:(N+m),
    nfeat = c(0,cum_nfeat),
    Method = rep("sample", N+m+1)
  )
  
  nfeat_sample <- cum_nfeat[N]
  
  bands <- data.frame(
    x = rep(N:(N+m), 2),
    means = c(nfeat_sample, ci_efpf$means + nfeat_sample,
              nfeat_sample, ci_mm$means + nfeat_sample),
    lbs = c(nfeat_sample, ci_efpf$lbs + nfeat_sample,
            nfeat_sample, ci_mm$lbs + nfeat_sample),
    ubs = c(nfeat_sample, ci_efpf$ubs + nfeat_sample,
            nfeat_sample, ci_mm$ubs + nfeat_sample),
    Method = c(rep("EFPF", m+1), rep("MM", m+1))
    
  )
  
  ggfig <- ggplot(data = bands, aes(x=x,y=means, group= Method) ) + 
    geom_line(aes(color = Method), linetype = "dashed") +
    geom_ribbon(aes(ymin = lbs, ymax = ubs, fill = Method), alpha = 0.1) +
    geom_line(data = obs_sample, aes(x, nfeat), color="black", linetype="solid", size=0.5) +
    geom_segment(aes(x = N, y = 0, xend = N, yend = max(ubs[m+1],cum_nfeat[N+m]) + 10), color="grey",
                 linetype="dashed", size=1) +
    xlab("# observations") + ylab("# distinct features") + theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0("N = ", N)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_x_continuous(breaks = pretty_breaks())
    
  
  return (ggfig)
}
 
  
#####################################################
#' Plot function for the observed trajectory 
#' 
#'
#' This function allows to plot the observed trajectory
#'
#' @param data_list [list] list of features in the whole sample
#'
#' @export
#' @import ggplot2
#' @import scales
#'
plot_trajectory <- function(data_list){
  
  L <- length(data_list)
  
  cum_nfeat <- sapply(1:L, function(n) length(unique(unlist(data_list[1:n]))))
  
  obs_sample <- data.frame(
    x = 0:L,
    nfeat = c(0,cum_nfeat)
  )
  
  
  ggfig <- ggplot(obs_sample, aes(x,nfeat) ) +
    geom_line(col = "black", linetype = "solid", size=0.5) +
    xlab("# observations") + ylab("# distinct features") + theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0("L = ", L)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_x_continuous(breaks = pretty_breaks())
  
  return (ggfig)
}

##############################################################

#' Plot Binary matrix of the feature observed for each individual
#'
#' This function allows to plot the features observed for each individual
#' following the order-of-appearance
#'
#' @param mat [matrix] binary matrix (#individuals x #features)
#' @param max_f [numeric] optional argument indicating the maximum
#' number of columns to plot
#' 
#' @export
#' 
#' @import ggplot2
#'
plot_binary_matrix <- function(mat, max_f = NULL){
  
  if(!is.null(max_f)){
    mat <- mat[,1:min(max_f, ncol(mat))]
  }
  
  mat <- mat[nrow(mat):1, ]
  
  colnames(mat) <- colnames(mat, do.NULL = FALSE, prefix = "f")
  rownames(mat) <- rownames(mat, do.NULL = FALSE, prefix = "ind")
  df <- reshape2::melt(mat)
  
  df$value <- factor(df$value)
  
  gg <- ggplot(df)
  
  # fill + legend, gray border
  gg <- gg + geom_tile(aes(x=Var2, y=Var1, fill=value), color="#7f7f7f")
  
  # custom fill colors
  gg <- gg + scale_fill_manual(values=c("white", "black"))
  
  # squares
  gg <- gg + coord_equal()
  
  # no labels
  gg <- gg + labs(x=NULL, y=NULL)
  
  # remove some chart junk
  gg <- gg + theme_bw()
  gg <- gg + theme(panel.grid=element_blank())
  gg <- gg + theme(panel.border=element_blank())
  gg <- gg + theme(axis.ticks =element_blank())
  gg <- gg + theme(axis.text.x = element_blank())
  gg <- gg + theme(axis.text.y = element_blank())
  gg <- gg + theme(legend.position="none")
  
  gg
  
}





