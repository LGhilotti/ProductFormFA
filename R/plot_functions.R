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
    lb = ci$lb,
    ub = ci$ub
  )
  
  ggplot(bands, aes(x, means)) + # ggplot2 plot with confidence intervals
    geom_line(col = "darkblue") +
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.1, fill = "darkblue") +
    xlab("m") + ylab(expression(K[m]^n)) + theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) +
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





