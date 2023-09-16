##
##
###########################################
## SubFunction of plot confidence band.
## Input is x, a vector of mean
## Input is LCL, a vector of lower bound
## Input is UCL, a vector of upper bound
## Output, plot a confidence band
##
##
conf.reg=function(x,LCL,UCL,...) polygon(c(x,rev(x)),c(LCL,rev(UCL)), ...)


##
##
###########################################
## Function of estimating individual-based "bootstrap community" in order to obtain estimated bootstrap s.e.
## Input is Spec, a vector of species abundances
## Output, a vector of estimated relative abundance
##
##
EstiBootComm.Ind <- function(Spec)
{
  Sobs <- sum(Spec > 0) 	#observed species
  n <- sum(Spec)		  	#sample size
  f1 <- sum(Spec == 1) 	#singleton 
  f2 <- sum(Spec == 2) 	#doubleton
  a <- ifelse(f1 == 0, 0, (n - 1) * f1 / ((n - 1) * f1 + 2 * f2) * f1 / n)
  b <- sum(Spec / n * (1 - Spec / n) ^ n)
  w <- a / b  			#adjusted factor for rare species in the sample
  f0.hat <- ceiling(ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2))	#estimation of unseen species via Chao1
  Prob.hat <- Spec / n * (1 - w * (1 - Spec / n) ^ n)					#estimation of relative abundance of observed species in the sample
  Prob.hat.Unse <- rep(2 * f2/((n - 1) * f1 + 2 * f2), f0.hat)		#estimation of relative abundance of unseen species in the sample
  return(c(Prob.hat, Prob.hat.Unse))									#Output: a vector of estimated relative abundance
}


##
##
###########################################
## Function of estimating sample-based "bootstrap community" in order to obtain estimated bootstrap s.e.
## Input is Spec, a vector of species abundances
## Input is T, number of samples
## Output, a vector of estimated detection probability
##
##
EstiBootComm.Sam <- function(Spec, T)
{
  Sobs <- sum(Spec > 0) 	#observed species
  Q1 <- sum(Spec == 1) 	#singleton 
  Q2 <- sum(Spec == 2) 	#doubleton
  a <- ifelse(Q1 == 0, 0, (T - 1) * Q1 / ((T - 1) * Q1 + 2 * Q2) * Q1 / T)
  b <- sum(Spec / T * (1 - Spec / T) ^ T)
  w <- a / b  			#adjusted factor for rare species in the sample
  Q0.hat <- ceiling(ifelse(Q2 == 0, (T - 1) / T * Q1 * (Q1 - 1) / 2, (T - 1) / T * Q1 ^ 2/ 2 / Q2))	#estimation of unseen species via Chao2
  Prob.hat <- Spec / T * (1 - w * (1 - Spec / T) ^ T)					#estimation of detection probability of observed species in the sample
  Prob.hat.Unse <- rep(2 * Q2/((T - 1) * Q1 + 2 * Q2), Q0.hat)		#estimation of detection probability of unseen species in the sample
  return(c(Prob.hat,  Prob.hat.Unse))									#Output: a vector of estimated detection probability
}


##
##
###########################################
## Estimation of interpolation and extrapolation of individual-based Hill number
## Input is x, a vector of species abundances
## Input is q, a numerical value of the order of Hill number
## Input is m, a integer vector of rarefaction/extrapolation sample size
## Output, a vector of estimated interpolation and extrapolation function of Hill number with order q
##
##
Dqhat.Ind <- function(x, q, m)
{
  x <- x[x > 0]
  n <- sum(x)
  
  fk.hat <- function(x, m){
    x <- x[x > 0]
    n <- sum(x)
    if(m <= n){
      Sub <- function(k)	sum(exp(lchoose(x, k) + lchoose(n - x, m -k) - lchoose(n, m)))
      sapply(1:m, Sub)
    }
    
    else {
      f1 <- sum(x == 1)
      f2 <- sum(x == 2)
      A <- ifelse(f2 > 0, (n-1)*f1/((n-1)*f1+2*f2), (n-1)*f1/((n-1)*f1+2))
      C.hat <- 1 - f1 / n * A
      p.hat <- x / n * C.hat			
      Sub <- function(k)	sum((choose(m, k) * p.hat^k * (1 - p.hat)^(m - k)) / (1 - (1 - p.hat)^n))
      sapply(1:m, Sub)
    }
  }
  
  D0.hat <- function(x, m){
    x <- x[x > 0]
    n <- sum(x)
    Sub <- function(m){
      if(m <= n){
        Fun <- function(x){
          if(x <= (n - m)) exp(lgamma(n - x + 1) + lgamma(n - m + 1) - lgamma(n - x - m + 1) - lgamma(n + 1))
          else 0
        }
        sum(1 - sapply(x, Fun))
      }
      else {
        Sobs <- sum(x > 0)
        f1 <- sum(x == 1)
        f2 <- sum(x == 2)
        f0.hat <- ifelse(f2 > 0, f1^2 /(2 * f2), f1 * (f1 - 1) / 2)
        ifelse(f1 ==0, Sobs ,Sobs + f0.hat * (1 - (1 - f1 / (n * f0.hat + f1)) ^ (m - n)))	
      }
    }
    sapply(m, Sub)
  }
  
  D1.hat <- function(x, m){
    x <- x[x > 0]
    n <- sum(x)
    Sub <- function(m){
      if(m < n){
        k <- 1:m
        exp(-sum(k / m * log(k / m) * fk.hat(x, m)))
      }
      else{
        #UE=sum(sapply(1:(n-1),function(k){sum(1/k*x/n*exp(lchoose(n-x,k)-lchoose(n-1,k)))}))
        UE <- sum(x/n*(digamma(n)-digamma(x)))
        f1 <- sum(x == 1)
        f2 <- sum(x == 2)
        A <- 1 - ifelse(f2 > 0, (n-1)*f1/((n-1)*f1+2*f2), (n-1)*f1/((n-1)*f1+2))
        #A=2*sum(x==2)/((n-1)*sum(x==1)+2*sum(x==2))
        B=sum(x==1)/n*(1-A)^(-n+1)*(-log(A)-sum(sapply(1:(n-1),function(k){1/k*(1-A)^k})))
        H.hat <- UE+B
        Hn.hat <- -sum(x / n * log(x / n))
        w <- (m - n) / m
        
        exp(w * H.hat + (1 - w) * Hn.hat)
        
      }
    }
    sapply(m, Sub)
  }
  
  D2.hat <- function(x, m){
    Sub <- function(m) 1 / (1 / m + (1 - 1 / m) * sum(x * (x - 1) / n / (n - 1)))
    sapply(m, Sub)
  }
  
  Dq.hat <- function(x, m){
    Sub <- function(m){
      k <- 1:m
      sum( (k / m)^q * fk.hat(x, m))^(1 / (1 - q))
    }
    sapply(m, Sub)
  }
  if(q == 0) D0.hat(x, m)
  else if(q == 1) D1.hat(x, m)
  else if(q == 2) D2.hat(x, m)
  else Dq.hat(x, m)
}


##
##
###########################################
## Estimation of interpolation and extrapolation of sample-based Hill number
## Input is y, a vector of species incidence-based frequency
## Input is T, number of samples
## Input is q, a numerical value of the order of Hill number
## Input is t, a integer vector of rarefaction/extrapolation sample size
## Output, a vector of estimated interpolation and extrapolation function of Hill number with order q
##
##
Dqhat.Sam <- function(y, T, q, t){
  y <- y[y > 0]
  U <- sum(y)
  
  Qk.hat <- function(y, T, t){
    y <- y[y > 0]
    U <- sum(y)
    if(t <= T){
      Sub <- function(k)	sum(exp(lchoose(y, k) + lchoose(T - y, t - k) - lchoose(T, t)))
      sapply(1:t, Sub)
    }
    
    else {
      Q1 <- sum(y==1)
      Q2 <- sum(y==2)
      A <- ifelse(Q2 > 0, (T-1)*Q1/((T-1)*Q1+2*Q2), (T-1)*Q1/((T-1)*Q1+2))
      C.hat <- 1 - Q1 / U * A
      p.hat <- y / T * C.hat			
      Sub <- function(k)	sum((choose(t, k) * p.hat^k * (1 - p.hat)^(t - k)) / (1 - (1 - p.hat)^T))
      sapply(1:t, Sub)
    }
  }
  
  D0.hat <- function(y, T, t){
    y <- y[y > 0]
    U <- sum(y)
    Sub <- function(t){
      if(t <= T){
        Fun <- function(y){
          if(y <= (T - t)) exp(lgamma(T - y + 1) + lgamma(T - t + 1) - lgamma(T - y - t + 1) - lgamma(T + 1))
          else 0
        }
        sum(1 - sapply(y, Fun))
      }
      else {
        Sobs <- sum(y > 0)
        Q1 <- sum(y==1)
        Q2 <- sum(y==2)
        A <- ifelse(Q2 > 0, (T-1)*Q1/((T-1)*Q1+2*Q2), (T-1)*Q1/((T-1)*Q1+2))
        C.hat <- 1 - Q1 / U * A
        Q0.hat <- ifelse(Q2 > 0, Q1^2 /(2 * Q2), Q1 * (Q1 - 1) / 2)
        ifelse(Q1 ==0, Sobs ,Sobs + Q0.hat * (1 - (1 - Q1 / (T * Q0.hat + Q1)) ^ (t - T)))	
      }
    }
    sapply(t, Sub)
  }
  
  D1.hat <- function(y, T, t){
    y <- y[y > 0]
    U <- sum(y)
    Sub <- function(t){
      U <- sum(y)
      if(t < T){
        k <- 1:t	
        Ut.hat <- t / T * U
        exp(-sum(k / Ut.hat * log(k / Ut.hat) * Qk.hat(y, T, t)))
      }
      else {
        UE <- sum(y/T*(digamma(T)-digamma(y)))
        Q1 <- sum(y == 1)
        Q2 <- sum(y == 2)
        A <- 1 - ifelse(Q2 > 0, (T-1)*Q1/((T-1)*Q1+2*Q2), (T-1)*Q1/((T-1)*Q1+2))
        B=sum(y==1)/T*(1-A)^(-T+1)*(-log(A)-sum(sapply(1:(T-1),function(k){1/k*(1-A)^k})))
        H.hat <- UE+B
        H.hat <- T/U*H.hat-log(T/U)
        
        Hn.hat <- -sum(y / U * log(y / U))
        w <- (t - T) / t
        exp(w * H.hat + (1 - w) * Hn.hat)
      }
    }
    sapply(t, Sub)
  }
  
  D2.hat <- function(y, T, t){
    Sub <- function(t) 1 / (1 / t * T / U + (1 - 1 / t) * sum(y * (y - 1) / U^2 / (1 - 1 / T)))
    sapply(t, Sub)
  }
  
  Dq.hat <- function(y, T, t){
    Sub <- function(t){
      k <- 1:t
      Ut.hat <- U * t / T
      sum( (k / Ut.hat)^q * Qk.hat(y, T, t))^(1 / (1 - q))
    }
    sapply(t, Sub)
  }
  if(q == 0) D0.hat(y, T, t)
  else if(q == 1) D1.hat(y, T, t)
  else if(q == 2) D2.hat(y, T, t)
  else Dq.hat(y, T, t)
}


##
##
###############################################
## Estimation of individual-based sample coverage function
## Input is x, a vector of species abundances
## Input is m, a integer vector of rarefaction/extrapolation sample size
## Output, a vector of estimated sample coverage function
##
##
Chat.Ind <- function(x, m)
{
  x <- x[x>0]
  n <- sum(x)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  A <- ifelse(f2 > 0, (n-1)*f1/((n-1)*f1+2*f2), (n-1)*f1/((n-1)*f1+2))
  Sub <- function(m){
    if(m < n) out <- 1-sum(x / n * exp(lchoose(n - x, m)-lchoose(n - 1, m)))
    if(m == n) out <- 1-f1/n*A
    if(m > n) out <- 1-f1/n*A^(m-n+1)
    out
  }
  sapply(m, Sub)		
}


##
##
###############################################
## Estimation of sample-based sample coverage function
## Input is y, a vector of species incidence-based frequency
## Input is T, number of samples
## Input is t, a integer vector of rarefaction/extrapolation sample size
## Output, a vector of estimated sample coverage function
##
##
Chat.Sam <- function(y, T, t){
  y <- y[y>0]
  U <- sum(y)
  Q1 <- sum(y == 1)
  Q2 <- sum(y == 2)
  A <- ifelse(Q2 > 0, (T-1)*Q1/((T-1)*Q1+2*Q2), (T-1)*Q1/((T-1)*Q1+2))
  Sub <- function(t){
    if(t < T) out <- 1 - sum(y / U * exp(lchoose(T - y, t) - lchoose(T - 1, t)))
    if(t == T) out <- 1 - Q1 / U * A
    if(t > T) out <- 1 - Q1 / U * A^(t - T + 1)
    out
  }
  sapply(t, Sub)		
}


##
##
###############################################
## Main program for individual-based data
## Input is Spec, a vector of species abundances
## Input is endpoint, a endpoint for extrapolation, default is double the reference sample size
## Input is Knots, a number of knots of computation, default is 40
## Input is rd, rounds the values to the specified number of decimal places, default is tens (-1) and other suggested argument is units (0) or hundreds (-2)
## Input is sd, calculate bootstrap standard error and 95% confidence interval; default is TRUE
## Input is nboot, the number of bootstrap resampling times, default is 200
## Output, a list of interpolation and extrapolation Hill number with order 0, 1, 2 and sample coverage 
##
##
iNEXT.Ind <- function(Spec, endpoint=2*sum(Spec), Knots=40, rd=-1, sd=TRUE, nboot=200)
{
  n <- sum(Spec)		  	#sample size
  m <- c(round(seq(1, sum(Spec)-0.1^rd, length=floor(Knots/2)-1), rd), sum(Spec), round(seq(sum(Spec)+0.1^rd, endpoint, length=floor(Knots/2)), rd))
  m <- c(1, m[-1])
  D0.hat <- Dqhat.Ind(Spec, q=0, m)
  D1.hat <- Dqhat.Ind(Spec, q=1, m)
  D2.hat <- Dqhat.Ind(Spec, q=2, m)
  Cov.hat <- Chat.Ind(Spec, m)
  if(sd==TRUE)
  {
    Prob.hat <- EstiBootComm.Ind(Spec)
    Abun.Mat <- rmultinom(nboot, n, Prob.hat)
    
    error.0 <-  qnorm(0.975) * apply(apply(Abun.Mat, 2, function(x) Dqhat.Ind(x, q=0, m)), 1, sd, na.rm=TRUE)
    left.0  <- D0.hat - error.0
    right.0 <- D0.hat + error.0
    
    error.1 <-  qnorm(0.975) * apply(apply(Abun.Mat, 2, function(x) Dqhat.Ind(x, q=1, m)), 1, sd, na.rm=TRUE)
    left.1  <- D1.hat - error.1
    right.1 <- D1.hat + error.1
    
    
    error.2 <-  qnorm(0.975) * apply(apply(Abun.Mat, 2, function(x) Dqhat.Ind(x, q=2, m)), 1, sd, na.rm=TRUE)
    left.2  <- D2.hat - error.2
    right.2 <- D2.hat + error.2
    
    error.C <-  qnorm(0.975) * apply(apply(Abun.Mat, 2, function(x) Chat.Ind(x, m)), 1, sd, na.rm=TRUE)
    left.C  <- Cov.hat - error.C
    right.C <- Cov.hat + error.C
    
    
    out.0 <- cbind("m"=m, "D0.hat"=D0.hat, "Norm.CI.Low"=left.0, "Norm.CI.High"=right.0, "Cov.hat"=Cov.hat)
    out.1 <- cbind("m"=m, "D1.hat"=D1.hat, "Norm.CI.Low"=left.1, "Norm.CI.High"=right.1, "Cov.hat"=Cov.hat)
    out.2 <- cbind("m"=m, "D2.hat"=D2.hat, "Norm.CI.Low"=left.2, "Norm.CI.High"=right.2, "Cov.hat"=Cov.hat)
    out.C <- cbind("m"=m, "Cov.hat"=Cov.hat, "Norm.CI.Low"=left.C, "Norm.CI.High"=right.C)
    out <- list("q=0"=out.0, "q=1"=out.1, "q=2"=out.2, "Cov"=out.C)
  }
  else
  {
    out.0 <- cbind("m"=m, "D0.hat"=D0.hat, "Cov.hat"=Cov.hat)
    out.1 <- cbind("m"=m, "D1.hat"=D1.hat, "Cov.hat"=Cov.hat)
    out.2 <- cbind("m"=m, "D2.hat"=D2.hat, "Cov.hat"=Cov.hat)
    out.C <- cbind("m"=m, "Cov.hat"=Cov.hat)
    out <- list("q=0"=out.0, "q=1"=out.1, "q=2"=out.2, "Cov"=out.C)
  }
  return(out)
}


##
##
###############################################
## Main program for individual-based data
## Input is Spec, a vector of species incidence-based frequency
## Input is T, number of samples
## Input is endpoint, a endpoint for extrapolation, default is double the reference sample size
## Input is Knots, a number of knots of computation, default is 40
## Input is rd, rounds the values to the specified number of decimal places, default is units (0) and other suggested argument is tens (-1) or hundreds (-2)
## Input is sd, calculate bootstrap standard error and 95% confidence interval; default is TRUE
## Input is nboot, the number of bootstrap resampling times, default is 200
## Output, a list of interpolation and extrapolation Hill number with order 0, 1, 2 and sample coverage 
##
##
iNEXT.Sam <- function(Spec, T, endpoint=2*T, Knots=40, rd=0, sd=TRUE, nboot=200)
{
  U <- sum(Spec)		  	#sample size
  t <- c(round(seq(1, T-0.1^rd, length=floor(Knots/2)-1), rd), T, round(seq(T+0.1^rd, endpoint, length=floor(Knots/2)), rd))
  t <- c(1, t[-1])
  D0.hat <- Dqhat.Sam(Spec, T, q=0, t)
  D1.hat <- Dqhat.Sam(Spec, T, q=1, t)
  D2.hat <- Dqhat.Sam(Spec, T, q=2, t)
  Cov.hat <- Chat.Sam(Spec, T, t)
  if(sd==TRUE)
  {
    Prob.hat <- EstiBootComm.Sam(Spec, T)
    Abun.Mat <- t(sapply(Prob.hat, function(p) rbinom(nboot, T, p)))
    
    error.0 <-  qnorm(0.975) * apply(apply(Abun.Mat, 2, function(y) Dqhat.Sam(y, T, q=0, t)), 1, sd, na.rm=TRUE)
    left.0  <- D0.hat - error.0
    right.0 <- D0.hat + error.0
    
    error.1 <-  qnorm(0.975) * apply(apply(Abun.Mat, 2, function(y) Dqhat.Sam(y, T, q=1, t)), 1, sd, na.rm=TRUE)
    left.1  <- D1.hat - error.1
    right.1 <- D1.hat + error.1
    
    
    error.2 <-  qnorm(0.975) * apply(apply(Abun.Mat, 2, function(y) Dqhat.Sam(y, T, q=2, t)), 1, sd, na.rm=TRUE)
    left.2  <- D2.hat - error.2
    right.2 <- D2.hat + error.2
    
    error.C <-  qnorm(0.975) * apply(apply(Abun.Mat, 2, function(y) Chat.Sam(y, T, t)), 1, sd, na.rm=TRUE)
    left.C  <- Cov.hat - error.C
    right.C <- Cov.hat + error.C
    
    
    out.0 <- cbind("t"=t, "D0.hat"=D0.hat, "Norm.CI.Low"=left.0, "Norm.CI.High"=right.0, "Cov.hat"=Cov.hat)
    out.1 <- cbind("t"=t, "D1.hat"=D1.hat, "Norm.CI.Low"=left.1, "Norm.CI.High"=right.1, "Cov.hat"=Cov.hat)
    out.2 <- cbind("t"=t, "D2.hat"=D2.hat, "Norm.CI.Low"=left.2, "Norm.CI.High"=right.2, "Cov.hat"=Cov.hat)
    out.C <- cbind("t"=t, "Cov.hat"=Cov.hat, "Norm.CI.Low"=left.C, "Norm.CI.High"=right.C)
    out <- list("q=0"=out.0, "q=1"=out.1, "q=2"=out.2, "Cov"=out.C)
  }
  else
  {
    out.0 <- cbind("t"=t, "D0.hat"=D0.hat, "Cov.hat"=Cov.hat)
    out.1 <- cbind("t"=t, "D1.hat"=D1.hat, "Cov.hat"=Cov.hat)
    out.2 <- cbind("t"=t, "D2.hat"=D2.hat, "Cov.hat"=Cov.hat)
    out.C <- cbind("t"=t, "Cov.hat"=Cov.hat)
    out <- list("q=0"=out.0, "q=1"=out.1, "q=2"=out.2, "Cov"=out.C)
  }
  return(out)
}

library(compiler)
iNEXT.Ind <- cmpfun(iNEXT.Ind)
iNEXT.Sam <- cmpfun(iNEXT.Sam)
##
##
###########################################
## Function of drawing figure
## input is out, a data frame of interpolation and extrapolation output. First column is the vector for x-axis, second column is for y-axis, and last two columns are pointwise confidence interval.
##
##

plot.iNEXT <- function(out, xlab=colnames(out)[1], ylab=colnames(out)[2], xlim=range(out[,1]), ylim=range(out[,-1]), col=1, ...)
{
  ref <- floor(nrow(out)/2)
  Inte <- as.data.frame(out[1:ref, ])
  Extr <- as.data.frame(out[ref:nrow(out), ])
  Mat <- rbind(Inte, Extr)
  
  if(ncol(Mat) < 4)
  {
    plot(0, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,...)
    lines(Inte, lty=1, lwd=2, col=col)
    lines(Extr, lty=2, lwd=2, col=col)
    points(Inte[ref,], pch=16, cex=1.5)		
  }
  
  else
  {
    conf.reg=function(x,LCL,UCL,...) polygon(c(x,rev(x)),c(LCL,rev(UCL)), ...)	#SubFunction of plot confidence band.
    plot(0, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab,...)
    conf.reg(Mat[,1], Mat$Norm.CI.Low, Mat$Norm.CI.High, col=adjustcolor(col, 0.25), border=NA)
    lines(Inte, lty=1, lwd=2, col=col)
    lines(Extr, lty=2, lwd=2, col=col)
    points(Inte[ref,], pch=16, cex=1.5, col=col)
  }
}

##
##
###########################################
## Example individual-based data, spiders abundance data collected by Sackett et al. (2011)
##
##
Girdled <- c(46, 22, 17, 15, 15, 9, 8, 6, 6, 4, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
Logged <- c(88, 22, 16, 15, 13, 10, 8, 8, 7, 7, 7, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

##
##
###########################################
## Example sample-based data, tropical ant species data collected by Longino and Colwell (2011)
## Note that first cell is number of total samples, and others are species incidence-based frequency.
##

## 50m
y50 <- c(599,rep(1,49),rep(2,23),rep(3,18),rep(4,14),rep(5,9),rep(6,10),rep(7,4),
         rep(8,8),rep(9,6),rep(10,2),rep(11,1),12,12,13,13,rep(14,5),15,15,
         rep(16,4),17,17,17,18,18,19,19,20,20,20,21,22,23,23,25,27,27,29,30,30,
         31,33,39,40,43,46,46,47,48,51,52,52,56,56,58,58,61,61,65,69,72,77,79,82,
         83,84,86,91,95,97,98,98,106,113,124,126,127,128,129,129,182,183,186,195,
         222,236,263,330)

##500m
y500 <- c(230,rep(1,71),rep(2,34),rep(3,12),rep(4,14),rep(5,9),rep(6,11),rep(7,8),
          rep(8,4),rep(9,7),rep(10,5),rep(11,2),12,12,12,13,13,13,13,14,14,15,
          16,16,17,17,17,17,18,19,20,21,21,23,24,25,25,25,26,27,30,31,31,32,32,
          33,34,36,37,38,38,38,38,39,39,41,42,43,44,45,46,47,49,52,52,53,54,56,
          60,60,65,73,78,123,131,133)

##1070m
y1070 <- c(150,rep(1,28),rep(2,16),rep(3,13),rep(4,3),rep(5,1),rep(6,3),rep(7,6),
           rep(8,1),rep(9,1),rep(10,1),rep(11,4),12,12,12,13,13,13,13,14,15,
           16,16,16,16,18,19,19,21,22,23,24,25,25,25,26,30,31,31,31,32,34,36,
           38,39,43,43,45,45,46,54,60,68,74,80,96,99)
##1500m
y1500 <- c(200,rep(1,13),rep(2,4),rep(3,2),rep(4,2),rep(5,4),rep(6,2),rep(9,4),
           rep(11,2),rep(17,2),18,19,23,23,24,25,25,25,29,30,32,33,43,50,53,
           73,74,76,79,113,144)

##2000m
y2000=c(200,1,2,2,3,4,8,8,13,15,19,23,34,59,80)