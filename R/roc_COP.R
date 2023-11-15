#' roc_COP
#'
#' Function to produce time-dependent ROC using Rotated 90 Clayton copula.
#'
#' @param max.x A value that specify what is the maximum X of our data.
#' @param time.t A vector of T that we are interested to produce its ROC.
#' @param mod A string that specify our joint model of biomarker X and time-to-event T.\cr
#'            mod = c('exp_exp','norm_exp','norm_weib')
#' @param params A list of estimated parameters.\cr
#'               'exp_exp': params = list( theta, lambda.t, mu.x).\cr
#'               'norm_exp': params = list(theta, lambda.t, mu.x, sigma.x).\cr
#'               'norm_weib': params = list(theta, scale.t, shape.t, mu.x, sigma.x).\cr
#'               theta: dependence measurement (tau (0 to 1);theta=2*tau/(1-tau))\cr
#'               lambda.t: scale parameter of T\cr
#'               scale.t: scale parameter of T (for weibull)\cr
#'               shape.t: shape parameter of T (for weibull)\cr
#'               mu.x: scale parameter of X\cr
#'               sigma.x: standard deviation of X
#' @param n.quant A number of how many quantile of Specificity that will be used in the ROC curve.
#' @export
#' @returns A dataframe of 4 columns (sensitivity, specificity, time, cutoff.x)
#' @examples
#' tau=0.5; theta=2*tau/(1-tau); n.quant=21; mu.x=1.4; lambda.t=0.2;
#' params=list(mu.x=mu.x,lambda.t=lambda.t,theta=theta)
#' roc_COP(max.x=max(df.COP$x), time.t=t, mod='exp_exp', params=params, n.quant=n.quant)



roc_COP <- function(max.x,time.t,mod,params, n.quant=10){

  # allocating memory & seed
  lambda.t = params$lambda.t; mu.x = params$mu.x; theta = params$theta
  sigma.x = params$sigma.x; shape.t = params$shape.t; scale.t = params$scale.t
  denom.se <- rep(NA,length(time.t)); denom.fpr <- rep(NA,length(time.t))
  ans <- data.frame('sensitivity'=rep(NA,n.quant*length(time.t)),
                    'specificity'=rep(seq(0,1,length.out=n.quant),length(time.t)),
                    'time'=sort(rep(time.t,n.quant)),
                    'cutoff.x'=rep(NA,n.quant*length(time.t)))

  set.seed(12345)

  if(mod=='exp_exp'){
    # calculating denominator of Sensitivity
    for (i in 1:length(time.t)){
      denom.se[i] <- pexp(time.t[i],rate=lambda.t)
    }

    # calculating denominator of Specificity
    ans$denom.se <- sort(rep(denom.se,n.quant))
    v <- ans$denom.se
    ans$denom.sp <- 1-ans$denom.se
    ans$prob.sp <- ans$specificity * ans$denom.sp + v

    # calculating inverse of Specificity (to find cutoff)
    u <- seq(0,1,length.out=10001)
    for (i in 1:length(denom.se)){
      u.check <- u + ((1-u)^-theta + denom.se[i]^-theta - 1)^(-1/theta)
      for (j in 1:length(ans$prob.sp)){
        if(denom.se[i] == ans$denom.se[j]){
          u.quant <- u[which.min(abs(u.check-ans$prob.sp[j]))]
          ans$cutoff.x[j] <- qexp(u.quant,rate=mu.x)
        }
      }
    }

    # calculating Sensitivity
    u <- pexp(ans$cutoff.x,rate=mu.x)
    ans$sensitivity <- ((1-u)^-theta + v^-theta - 1)^-(1/theta)/v
  }

  else if(mod=='norm_exp'){
    # calculating denominator of Sensitivity
    for (i in 1:length(time.t)){
      denom.se[i] <- pexp(time.t[i],rate=lambda.t)
    }

    # calculating denominator of Specificity
    ans$denom.se <- sort(rep(denom.se,n.quant))
    v <- ans$denom.se
    ans$denom.sp <- 1-ans$denom.se
    ans$prob.sp <- ans$specificity * ans$denom.sp + v

    # calculating inverse of Specificity (to find cutoff)
    u <- seq(0,1,length.out=10001)
    for (i in 1:length(denom.se)){
      u.check <- u + ((1-u)^-theta + denom.se[i]^-theta - 1)^(-1/theta)
      for (j in 1:length(ans$prob.sp)){
        if(denom.se[i] == ans$denom.se[j]){
          u.quant <- u[which.min(abs(u.check-ans$prob.sp[j]))]
          ans$cutoff.x[j] <- qnorm(u.quant,mean=mu.x, sd=sigma.x)
        }
      }
    }

    # calculating Sensitivity
    u <- pnorm(ans$cutoff.x,mean=mu.x, sd=sigma.x)
    ans$sensitivity <- ((1-u)^-theta + v^-theta - 1)^-(1/theta)/v
  }

  else if(mod=='norm_weib'){
    # calculating denominator of Sensitivity
    for (i in 1:length(time.t)){
      denom.se[i] <- pweibull(time.t[i],scale=1/scale.t,shape=shape.t)
    }

    # calculating denominator of Specificity
    ans$denom.se <- sort(rep(denom.se,n.quant))
    v <- ans$denom.se
    ans$denom.sp <- 1-ans$denom.se
    ans$prob.sp <- ans$specificity * ans$denom.sp + v

    # calculating inverse of Specificity (to find cutoff)
    u <- seq(0,1,length.out=10001)
    for (i in 1:length(denom.se)){
      u.check <- u + ((1-u)^-theta + denom.se[i]^-theta - 1)^(-1/theta)
      for (j in 1:length(ans$prob.sp)){
        if(denom.se[i] == ans$denom.se[j]){
          u.quant <- u[which.min(abs(u.check-ans$prob.sp[j]))]
          ans$cutoff.x[j] <- qnorm(u.quant,mean=mu.x, sd=sigma.x)
        }
      }
    }

    # calculating Sensitivity
    u <- pnorm(ans$cutoff.x,mean=mu.x, sd=sigma.x)
    ans$sensitivity <- ((1-u)^-theta + v^-theta - 1)^-(1/theta)/v
  }
  ans$sensitivity[which(ans$sensitivity>1)] <- 1
  for(i in 1:length(time.t)){
    ans[n.quant*(i-1)+1,] <- c(1,0,time.t[i],NA,NA,NA,NA)
    ans[n.quant*i,] <- c(0,1,time.t[i],NA,NA,NA,NA)
  }
  ans <- ans[,c(1,2,3,4)]
  return(ans)
}
