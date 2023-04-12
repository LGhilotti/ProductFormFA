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
  
  
  ggplot(bands, aes(x,means) ) +
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
  
  
  ggplot(bands, aes(x,means) ) +
    geom_line(col = "red", linetype = "dashed") +
    geom_ribbon(aes(ymin = lbs, ymax = ubs), alpha = 0.1, fill = "red") +
    geom_line(data = obs_sample, aes(x, nfeat), color="black", linetype="solid", size=0.5) +
    geom_segment(aes(x = N, y = 0, xend = N, yend = ubs[m+1] + 10), color="grey",
                 linetype="dashed", size=1) +
    xlab("# observations") + ylab("# distinct features") + theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste0("N = ", N)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_x_continuous(breaks = pretty_breaks())
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





