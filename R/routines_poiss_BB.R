library(Rcpp)
library(ggplot2)
library(scales)

sourceCpp("src/Poiss_BB.cpp")

plot_Kmn_poiss_BB <- function(alpha, theta, m, n, lambda, lev) {
  means <- mean_kmn_all_poiss_BB(alpha, theta, m, n, lambda)
  ub <- qpois(lev + (1 - lev) / 2, means)
  lb <- qpois((1 - lev) / 2, means)
  bands <- data.frame(
    x = 1:m,
    means = means,
    lb = lb,
    ub = ub
  )

  ggplot(bands, aes(x, means)) + # ggplot2 plot with confidence intervals
    geom_line(col = "darkblue") +
    geom_ribbon(aes(ymin = lb, ymax = ub), alpha = 0.1, fill = "darkblue") +
    xlab("m") + ylab(expression(K[m]^n)) + theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_x_continuous(breaks = pretty_breaks())
  # scale_y_continuous(breaks = seq(0, max(ub)+1, by = round(max(ub)/10))) +
}
