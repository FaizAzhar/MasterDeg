#' roc_multisn
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
#' roc_multisn(max.x=max(df.COP$x), time.t=t, mod='exp_exp', params=params, n.quant=n.quant)
#'



roc_multisn <- function(max.x,time.t,mod,params,n.sample=100000, cutoff.n=10001, n.quant=10, low.lim=NULL, upp.lim=NULL){

  # allocating memory & seed
  mu.sn = params$mu.sn; sigma.sn = params$sigma.sn; skew.sn = params$skew.sn
  ans <- data.frame('sensitivity'=rep(NA,n.quant*length(time.t)),
                        'specificity'=rep(NA,n.quant*length(time.t)),
                        'time'=sort(rep(time.t,n.quant)),
                        'cutoff.x'=rep(seq(0,max.x,length.out=n.quant),length(time.t)))
  denom.se <- NULL
  numer.sp <- NULL
  numer.se <- NULL
  if(mod=="multi_snorm"){
    # calculating denominator of Sensitivity
    if(mod=="multi_snorm"){
      f.joint <- function(vars){
        x <- vars[1]
        t.T <- vars[2]
        dmsn(x=c(x,t.T), xi=mu.sn, Omega=sigma.sn, alpha=skew.sn)
      }
    }

    for (i in 1:length(time.t)){
      double.integ <- hcubature(f.joint, lowerLimit=c(low.lim[1],low.lim[2]), upperLimit=c(upp.lim[1],time.t[i]))$integral
      res <- rep(double.integ, n.quant)
      denom.se <- append(denom.se, res)
    }

    # calculating denominator of Specificity
    ans$denom.se <- denom.se
    ans$denom.sp <- 1-ans$denom.se

    # calculating Specificity's numerator
    for (i in 1:nrow(ans)){
      double.integ <- hcubature(f.joint, c(low.lim[1],ans$time[i]), c(ans$cutoff.x[i],upp.lim[2]))$integral
      numer.sp <- append(numer.sp, double.integ)
      print(paste0("Calculating Specificity's numerator: iter = ",i))
    }
    ans$numer.sp <- numer.sp

    # calculating Sensitivity's numerator
    for (i in 1:nrow(ans)){
      double.integ <- hcubature(f.joint, c(ans$cutoff.x[i],low.lim[2]), c(upp.lim[1],ans$time[i]))$integral
      numer.se <- append(numer.se, double.integ)
      print(paste0("Calculating Sensitivity's numerator: iter = ",i))
    }
    ans$numer.se <- numer.se

    ans$sensitivity <- ans$numer.se / ans$denom.se
    ans$specificity <- ans$numer.sp / ans$denom.sp
  }

  ans$sensitivity[which(ans$sensitivity>1)] <- 1
  for(i in 1:length(time.t)){
    ans <- add_row(ans,sensitivity=1,specificity=0,time=time.t[i],cutoff.x=NA,denom.se=NA,denom.sp=NA,numer.sp=NA,numer.se=NA,.before=(n.quant*(i-1)+2*i-1))
    ans <- add_row(ans,sensitivity=0,specificity=1,time=time.t[i],cutoff.x=NA,denom.se=NA,denom.sp=NA,numer.sp=NA,numer.se=NA,.before=(n.quant*(i)+2*i))
  }

  ans <- ans[,c(1,2,3,4)]
  return(ans)
}
