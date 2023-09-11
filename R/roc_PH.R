roc_PH <- function(max.x,time.t,mod,params,n.sample=100000, cutoff.n=10001, n.quant=10){

  # allocating memory & seed
  lambda.t = params$lambda.t; mu.x = params$mu.x; beta = params$beta;
  sigma.x = params$sigma.x; shape.t = params$shape.t; scale.t = params$scale.t
  denom.se <- rep(NA,length(time.t)); denom.fpr <- rep(NA,length(time.t))
  ans <- data.frame('sensitivity'=rep(NA,n.quant*length(time.t)),
                    'specificity'=rep(seq(0,1,length.out=n.quant),length(time.t)),
                    'time'=sort(rep(time.t,n.quant)),
                    'cutoff.x'=rep(NA,n.quant*length(time.t)))
  x <- rep(NA, n.sample); t <- rep(NA, n.sample);
  pb.cdf <- txtProgressBar(min=0, max=cutoff.n, style=3,
                           width=50, char="*")
  pb.overall <- txtProgressBar(min=0,max=n.quant*length(time.t),
                               style=3, width=50, char="*")
  set.seed(12345)

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
    cat("\nCreating CDF table: \n")
    for (t in time.t){
      for(i in 1: cutoff.n){
        res[i] <- cdf_trapz(f.c,b=cutoff[i])
        if(is.na(res[i])) res[i] <- 0
        setTxtProgressBar(pb.cdf,i)

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
        setTxtProgressBar(pb.overall,i)
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
      setTxtProgressBar(pb.overall,i)
    }
    close(pb.overall)
    ans$sensitivity <- ans$sensitivity/ans$denom.se
    print(Sys.time()-start.time)
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
    cat("\nCreating CDF table: \n")
    start.time <- Sys.time()
    for (t in time.t){
      for(i in 1: cutoff.n){
        res[i] <- cdf_trapz(f.c,b=cutoff[i])
        if(is.na(res[i])) res[i] <- 0
        setTxtProgressBar(pb.cdf,i)
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
        setTxtProgressBar(pb.overall,i)
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
      setTxtProgressBar(pb.overall,i)
    }
    close(pb.overall)
    ans$sensitivity <- ans$sensitivity/ans$denom.se
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
    cat("\nCreating CDF table: \n")
    start.time <- Sys.time()
    for (t in time.t){
      for(i in 1: cutoff.n){
        res[i] <- cdf_trapz(f.c,b=cutoff[i])
        if(is.na(res[i])) res[i] <- 0
        setTxtProgressBar(pb.cdf,i)
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
        setTxtProgressBar(pb.overall,i)
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
      setTxtProgressBar(pb.overall,i)
    }
    close(pb.overall)
    ans$sensitivity <- ans$sensitivity/ans$denom.se
    print(Sys.time()-start.time)
  }
  ans$sensitivity[which(ans$sensitivity>1)] <- 1
  for(i in 1:length(time.t)){
    ans[n.quant*(i-1)+1,] <- c(1,0,time.t[i],NA,NA,NA,NA)
    ans[n.quant*i,] <- c(0,1,time.t[i],NA,NA,NA,NA)
  }
  ans <- ans[,c(1,2,3,4)]
  return(ans)
}
