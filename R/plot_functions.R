#' Plot function for the credible intervals of Kmn
#'
#' This function allows to plot the credible intervals of Kmn
#'
#' @param ci [list] it contains means, upper-bounds and lower-bounds of the credible intervals
#'
#' @export
#' @import ggplot2, scales
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
