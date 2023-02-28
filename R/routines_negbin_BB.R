library(Rcpp)
library(ggplot2)
library(scales)

sourceCpp("src/NegBin_BB.cpp")


plot_Kmn_negbin_BB <- function(alpha,theta, m, n, Kn, nstar, p, lev){
  
  pbars <- p_kmn_all_negbin_BB(alpha,theta, m, n, p)
  
  means <- (nstar + Kn)*(1-pbars)/pbars 
  
  ub <- qnbinom(lev+(1-lev)/2, nstar+Kn, pbars)
  
  lb <- qnbinom((1-lev)/2, nstar+Kn, pbars)
  
  bands <- data.frame(x=1:m,
                      means = means,
                      lb = lb,
                      ub = ub)
  
  ggplot(bands, aes(x, means)) +        # ggplot2 plot with confidence intervals
    geom_point() +
    geom_errorbar(aes(ymin = lb, ymax = ub)) +
    labs(title=paste0("CI of Kmn for BB with NB(",nstar,",",p, "), with n=",n,", Kn=",Kn),
         x ="m", y = "Kmn") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(breaks = pretty_breaks()) +
    scale_x_continuous(breaks= pretty_breaks())

}
