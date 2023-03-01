library(Rcpp)
library(ggplot2)
library(scales)

sourceCpp("src/Gamma_IBP.cpp")


plot_Kmn_gamma_IBP <- function(alpha, theta, m, n, Kn, a, b, lev) {
  pbars <- p_kmn_all_gamma_IBP(alpha, theta, m, n, b)
  means <- (a + Kn) * (1 - pbars) / pbars
  ub <- qnbinom(lev + (1 - lev) / 2, a + Kn, pbars)
  lb <- qnbinom((1 - lev) / 2, a + Kn, pbars)
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
}
