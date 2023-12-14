#' roc_PH
#'
#' Function to produce time-dependent ROC using Proportional Hazard Model.
#'
#' @param max.x A value that specify what is the maximum X of our data.
#' @param time.t A vector of T that we are interested to produce its ROC.
#' @param mod A string that specify our joint model of biomarker X and time-to-event T. \cr
#'            mod = c('exp_exp','norm_exp','norm_weib')
#' @param params A list of estimated parameters.\cr
#'               'exp_exp': params = list(beta, lambda.t, mu.x).\cr
#'               'norm_exp': params = list(beta, lambda.t, mu.x, sigma.x).\cr
#'               'norm_weib': params = list(beta, scale.t, shape.t, mu.x, sigma.x).\cr
#'               beta: association measurement\cr
#'               lambda.t: scale parameter of T\cr
#'               scale.t: scale parameter of T (for weibull)\cr
#'               shape.t: shape parameter of T (for weibull)\cr
#'               mu.x: scale parameter of X\cr
#'               sigma.x: standard deviation of X
#' @param n.sample A value of sample size that will be used in Monte Carlo Integration method (to find inverse of Specificity).
#' @param cutoff.n A value of how many biomarker subdivision is used (to find inverse of Specificity).
#' @param n.quant A number of how many quantile of Specificity that will be used in the ROC curve.
#' @export
#' @returns A dataframe of 4 columns (sensitivity, specificity, time, cutoff.x)
#' @examples
#' beta=1.2; n.quant=21; mu.x=1.4; lambda.t=0.2;
#' params=list(mu.x=mu.x,lambda.t=lambda.t,beta=beta)
#' roc_PH(max.x=max(df.COP$x), time.t=t, mod='exp_exp', params=params, n.quant=n.quant)
#'



roc_PH <- function(max.x, min.x=NULL, time.t,mod,params,n.sample=100000, cutoff.n=10001, n.quant=10, low.lim=NULL, upp.lim=NULL){

  # allocating memory & seed
  lambda.t = params$lambda.t; mu.x = params$mu.x; beta = params$beta;
  sigma.x = params$sigma.x; shape.t = params$shape.t; scale.t = params$scale.t;
  skew.x = params$skew.x; mu.t = params$mu.t; sigma.t = params$sigma.t; skew.t = params$skew.t
  mu.sn = params$mu.sn; sigma.sn = params$sigma.sn; skew.sn = params$skew.sn
  denom.se <- rep(NA,length(time.t)); denom.fpr <- rep(NA,length(time.t))

  if(mod!="snorm_snorm"){
    ans <- data.frame('sensitivity'=rep(NA,n.quant*length(time.t)),
                    'specificity'=rep(seq(0,1,length.out=n.quant),length(time.t)),
                    'time'=sort(rep(time.t,n.quant)),
                    'cutoff.x'=rep(NA,n.quant*length(time.t)))
  x <- rep(NA, n.sample); t <- rep(NA, n.sample);
  cat("Creating CDF table: \n")
  pb.cdf <- txtProgressBar(min=0,max=1,
                           style=3, width=50, char="*")
  pb.overall <- txtProgressBar(min=0,max=1,
                               style=3, width=50, char="*")
  set.seed(12345)
  }

  if(mod=='exp_exp'){
    # calculating denominator of Sensitivity
    u <- runif(n.sample,0,0.99)
    x <- u/(1-u)
    denom <- rep(NA,n.sample)
    for (i in 1:length(time.t)){
      t <- runif(n.sample,0,time.t[i])
      denom <- lambda.t*mu.x*exp(x*beta)*exp(-exp(beta*x)*lambda.t*t)*exp(-mu.x*x)/(1-u)^2
      denom.se[i] <- time.t[i]*(0.99)/n.sample*sum(denom)
    }

    # calculating denominator of Specificity
    ans$denom.se <- sort(rep(denom.se,n.quant))
    ans$denom.sp <- 1-ans$denom.se
    ans$prob.sp <- ans$specificity * ans$denom.sp

    # calculating inverse of Specificity (to find cutoff) using trapezoidal rule
    f.c <- function(c) mu.x*exp(-exp(beta*c)*lambda.t*t)*exp(-mu.x*c)
    # preparing CDF table splitted by time
    cutoff <- seq(0,max.x,length.out=cutoff.n)
    res <- rep(NA,cutoff.n)
    cdf.df <- data.frame('C'=cutoff)

    start.time <- Sys.time()
    for (j in 1:length(time.t)){
      t <- time.t[j]
      for(i in 1: cutoff.n){
        res[i] <- cdf_trapz(f.c,b=cutoff[i])
        if(is.na(res[i])) res[i] <- 0
        setTxtProgressBar(pb.cdf,(cutoff.n*(j-1)+i)/(cutoff.n*length(time.t)))
      }
      cdf.df[, ncol(cdf.df)+1] <- res
    }
    close(pb.cdf)
    print(Sys.time()-start.time)

    # inverse CDF
    start.time <- Sys.time()
    cat("\nFinding inverse CDF: \n")
    for(t in 1:length(time.t)){
      for(i in 1:n.quant){
        ans$cutoff.x[(t-1)*n.quant+i] <- cdf.df$C[which.min(abs(ans$prob.sp[(t-1)*n.quant+i]-cdf.df[,t+1]))]
        setTxtProgressBar(pb.overall,(i+(t-1)*n.quant)/(n.quant*length(time.t)))
      }
    }
    close(pb.overall)
    print(Sys.time()-start.time)

    # calculating Sensitivity
    start.time <- Sys.time()
    x <- rexp(n.sample,rate=mu.x)
    cat("\nCalculating Sensitivity: \n")
    indicator.f <- function(x,c,t){if(x>c) return(1-exp(-exp(x*beta)*lambda.t*t)) else return(0)}
    for (i in 1:length(ans$cutoff.x)){
      ans$sensitivity[i] <- mean(unlist(sapply(x,indicator.f, c=ans$cutoff.x[i], t=ans$time[i])))
      setTxtProgressBar(pb.overall,i/length(ans$cutoff.x))
    }
    close(pb.overall)
    cat("\n")
    print(Sys.time()-start.time)
    ans$sensitivity <- ans$sensitivity/ans$denom.se

  }

  else if(mod=='norm_exp'){
    # calculating denominator of Sensitivity
    u <- runif(n.sample,-0.99,0.99)
    x <- u/(1-u^2)
    denom <- rep(NA,n.sample)
    for (i in 1:length(time.t)){
      t <- runif(n.sample,0,time.t[i])
      denom <- lambda.t*exp(x*beta)*exp(-exp(beta*x)*lambda.t*t)*exp(-(x-mu.x)^2/2/sigma.x^2)/sqrt(2*pi*sigma.x^2)/(1-u^2)^2*(1+u^2)
      denom.se[i] <- time.t[i]*(2*0.99)/n.sample*sum(denom)
    }

    # calculating denominator of Specificity
    ans$denom.se <- sort(rep(denom.se,n.quant))
    ans$denom.sp <- 1-ans$denom.se
    ans$prob.sp <- ans$specificity * ans$denom.sp

    # calculating inverse of Specificity (to find cutoff) using trapezoidal rule
    f.c <- function(c) (1/sqrt(2*pi*sigma.x^2))*exp(-(c-mu.x)^2/(2*sigma.x^2))*exp(-exp(c*beta)*lambda.t*t)
    # preparing CDF table splitted by time
    cutoff <- seq(0,max.x,length.out=cutoff.n)
    res <- rep(NA,cutoff.n)
    cdf.df <- data.frame('C'=cutoff)
    start.time <- Sys.time()
    for (j in 1:length(time.t)){
      t <- time.t[j]
      for(i in 1: cutoff.n){
        res[i] <- cdf_trapz(f.c,b=cutoff[i])
        if(is.na(res[i])) res[i] <- 0
        setTxtProgressBar(pb.cdf,(cutoff.n*(j-1)+i)/(cutoff.n*length(time.t)))
      }
      cdf.df[, ncol(cdf.df)+1] <- res
    }
    close(pb.cdf)
    print(Sys.time()-start.time)

    # inverse CDF
    start.time <- Sys.time()
    cat("Finding inverse CDF: \n")
    for(t in 1:length(time.t)){
      for(i in 1:n.quant){
        ans$cutoff.x[(t-1)*n.quant+i] <- cdf.df$C[which.min(abs(ans$prob.sp[(t-1)*n.quant+i]-cdf.df[,t+1]))]
        setTxtProgressBar(pb.overall,(i+(t-1)*n.quant)/(n.quant*length(time.t)))
      }
    }
    close(pb.overall)
    print(Sys.time()-start.time)

    # calculating Sensitivity
    x <- rnorm(n.sample,mean=mu.x,sd=sigma.x)
    start.time <- Sys.time()
    cat("Calculating Sensitivity: \n")
    indicator.f <- function(x,c,t){if(x>c) return(1-exp(-exp(x*beta)*lambda.t*t)) else return(0)}
    for (i in 1:length(ans$cutoff.x)){
      ans$sensitivity[i] <- mean(unlist(sapply(x,indicator.f, c=ans$cutoff.x[i], t=ans$time[i])))
      setTxtProgressBar(pb.overall,i/length(ans$cutoff.x))
    }
    close(pb.overall)
    ans$sensitivity <- ans$sensitivity/ans$denom.se
    cat("\n")
    print(Sys.time()-start.time)
  }

  else if(mod=='norm_weib'){
    # calculating denominator of Sensitivity
    u <- runif(n.sample,-0.99,0.99)
    x <- u/(1-u^2)
    denom <- rep(NA,n.sample)
    for (i in 1:length(time.t)){
      t <- runif(n.sample,0,time.t[i])
      denom <- exp(x*beta)*scale.t*shape.t*(scale.t*t)^(shape.t-1)*exp(-exp(beta*x)*(scale.t*t)^shape.t)*exp(-(x-mu.x)^2/2/sigma.x^2)/sqrt(2*pi*sigma.x^2)/(1-u^2)^2*(1+u^2)
      denom.se[i] <- time.t[i]*(2*0.99)/n.sample*sum(denom)
    }

    # calculating denominator of Specificity
    ans$denom.se <- sort(rep(denom.se,n.quant))
    ans$denom.sp <- 1-ans$denom.se
    ans$prob.sp <- ans$specificity * ans$denom.sp

    # calculating inverse of Specificity (to find cutoff) using trapezoidal rule
    f.c <- function(c) (1/sqrt(2*pi*sigma.x^2))*exp(-(c-mu.x)^2/(2*sigma.x^2))*exp(-exp(c*beta)*(scale.t*t)^shape.t)
    # preparing CDF table splitted by time
    cutoff <- seq(0,max.x,length.out=cutoff.n)
    res <- rep(NA,cutoff.n)
    cdf.df <- data.frame('C'=cutoff)
    start.time <- Sys.time()
    for (j in 1:length(time.t)){
      t <- time.t[j]
      for(i in 1: cutoff.n){
        res[i] <- cdf_trapz(f.c,b=cutoff[i])
        if(is.na(res[i])) res[i] <- 0
        setTxtProgressBar(pb.cdf,(cutoff.n*(j-1)+i)/(cutoff.n*length(time.t)))
      }
      cdf.df[, ncol(cdf.df)+1] <- res
    }
    close(pb.cdf)
    print(Sys.time()-start.time)

    # inverse CDF
    start.time <- Sys.time()
    cat("\nFinding inverse CDF: \n")
    for(t in 1:length(time.t)){
      for(i in 1:n.quant){
        ans$cutoff.x[(t-1)*n.quant+i] <- cdf.df$C[which.min(abs(ans$prob.sp[(t-1)*n.quant+i]-cdf.df[,t+1]))]
        setTxtProgressBar(pb.overall,(i+(t-1)*n.quant)/(n.quant*length(time.t)))
      }
    }
    close(pb.overall)
    print(Sys.time()-start.time)

    # calculating Sensitivity
    x <- rnorm(n.sample,mean=mu.x,sd=sigma.x)
    start.time <- Sys.time()
    cat("\nCalculating Sensitivity: \n")
    indicator.f <- function(x,c,t){if(x>c) return(1-exp(-exp(x*beta)*(scale.t*t)^shape.t)) else return(0)}
    for (i in 1:length(ans$cutoff.x)){
      ans$sensitivity[i] <- mean(unlist(sapply(x,indicator.f, c=ans$cutoff.x[i], t=ans$time[i])))
      setTxtProgressBar(pb.overall,i/length(ans$cutoff.x))
    }
    close(pb.overall)
    ans$sensitivity <- ans$sensitivity/ans$denom.se
    cat("\n")
    print(Sys.time()-start.time)
  }
  else if(mod=='snorm_snorm'){
    # calculating denominator of Sensitivity
    ans <- data.frame('sensitivity'=rep(NA,n.quant*length(time.t)),
                      'specificity'=rep(NA,n.quant*length(time.t)),
                      'time'=sort(rep(time.t,n.quant)),
                      'cutoff.x'=rep(seq(min.x,max.x,length.out=n.quant),length(time.t)))

    denom.se <- NULL
    f.joint <- function(vars){
    x <- vars[2]
    t.T <- vars[1]
    exp(x*beta)*((1-psn(t.T, xi=mu.t, omega=sigma.t, alpha=skew.t))^(exp(x*beta)-1))*dsn(t.T, xi=mu.t, omega=sigma.t, alpha=skew.t)*dsn(x, xi=mu.x, omega=sigma.x, alpha=skew.x)
    }

    for (i in 1:length(time.t)){
      double.integ <- hcubature(f.joint, lowerLimit=c(low.lim[1],low.lim[2]), upperLimit=c(time.t[i],upp.lim[2]))$integral
      res <- rep(double.integ, n.quant)
      denom.se <- append(denom.se, res)
    }

    # calculating denominator of Specificity
    ans$denom.se <- denom.se
    ans$denom.sp <- 1-ans$denom.se

    # calculating Specificity's numerator
    numer.sp <- NULL
    for (i in 1:nrow(ans)){
      double.integ <- hcubature(f.joint, c(ans$time[i],low.lim[2]), c(upp.lim[1],ans$cutoff.x[i]))$integral
      numer.sp <- append(numer.sp, double.integ)
      print(paste0("Calculating Specificity's numerator: iter = ",i))
    }
    ans$numer.sp <- numer.sp

    # calculating Sensitivity's numerator
    numer.se <- NULL
    for (i in 1:nrow(ans)){
      double.integ <- hcubature(f.joint, c(low.lim[1],ans$cutoff.x[i]), c(ans$time[i],upp.lim[2]))$integral
      numer.se <- append(numer.se, double.integ)
      print(paste0("Calculating Sensitivity's numerator: iter = ",i))
    }
    ans$numer.se <- numer.se

    ans$sensitivity <- ans$numer.se / ans$denom.se
    ans$specificity <- ans$numer.sp / ans$denom.sp
  }
  if(mod!="snorm_snorm"){
    ans$sensitivity[which(ans$sensitivity>1)] <- 1
    for(i in 1:length(time.t)){
      ans[n.quant*(i-1)+1,] <- c(1,0,time.t[i],NA,NA,NA,NA)
      ans[n.quant*i,] <- c(0,1,time.t[i],NA,NA,NA,NA)
    }
  }else{
    ans$sensitivity[which(ans$sensitivity>1)] <- 1
    for(i in 1:length(time.t)){
      ans <- add_row(ans,sensitivity=1,specificity=0,time=time.t[i],cutoff.x=NA,denom.se=NA,denom.sp=NA,numer.sp=NA,numer.se=NA,.before=(n.quant*(i-1)+2*i-1))
      ans <- add_row(ans,sensitivity=0,specificity=1,time=time.t[i],cutoff.x=NA,denom.se=NA,denom.sp=NA,numer.sp=NA,numer.se=NA,.before=(n.quant*(i)+2*i))
    }
  }
  # ans <- ans[,c(1,2,3,4)]
  return(ans)
}
